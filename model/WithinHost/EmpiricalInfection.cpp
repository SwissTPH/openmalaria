/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
 OpenMalaria is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or (at
 your option) any later version.
 
 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include "WithinHost/Empirical.h"
#include "util/gsl.h"
#include "util/errors.hpp"
#include "util/CommandLine.hpp"
#include "util/ModelOptions.hpp"

#include <sstream>
#include <fstream>
#include <string>
#include <cmath>

namespace OM { namespace WithinHost {
    
// -----  static class members (variables & functions)  -----

double EmpiricalInfection::_maximumPermittedAmplificationPerCycle;
double EmpiricalInfection::_subPatentLimit;
double EmpiricalInfection::_alpha1;
double EmpiricalInfection::_alpha2;	
double EmpiricalInfection::_alpha3;
double EmpiricalInfection::_mu1;	
double EmpiricalInfection::_mu2;	
double EmpiricalInfection::_mu3;
double EmpiricalInfection::_sigma0_res;	
double EmpiricalInfection::_sigmat_res;
double EmpiricalInfection::_mu_beta1[_maximumDurationInDays];
double EmpiricalInfection::_sigma_beta1[_maximumDurationInDays];
double EmpiricalInfection::_mu_beta2[_maximumDurationInDays];
double EmpiricalInfection::_sigma_beta2[_maximumDurationInDays];
double EmpiricalInfection::_mu_beta3[_maximumDurationInDays];
double EmpiricalInfection::_sigma_beta3[_maximumDurationInDays];
double EmpiricalInfection::_inflationMean;
double EmpiricalInfection::_inflationVariance;
double EmpiricalInfection::_extinctionLevel;
double EmpiricalInfection::_overallMultiplier;


void EmpiricalInfection::initParameters(){
  if (Global::interval != 1)
    throw util::xml_scenario_error ("EmpiricalInfection only supports using an interval of 1");
  
  // alpha1 corresponds to 1 day before first patent, alpha2 2 days before first patent etc.
  _alpha1=0.2647;
  _alpha2=2.976;
  _alpha3=0.9181;
  _mu1=6.08e-04;
  _mu2=0.624;
  _mu3=0.3064;
  _sigma0_res=0.9998;
  _sigmat_res=0.002528;
/* The following variables are assigned separately for each infection to enable optimisation of their values
*/    
  _inflationMean= 1.09635;
  _inflationVariance= 0.172029;
  _extinctionLevel= 0.0100976;
  _overallMultiplier= 0.697581;
  _subPatentLimit=10.0/_overallMultiplier; 
  _maximumPermittedAmplificationPerCycle=1000.0;
  fstream f_autoRegressionParameters(util::CommandLine::lookupResource("autoRegressionParameters.csv").c_str(),ios::in);
  if (!f_autoRegressionParameters.is_open())
    throw runtime_error ("file not found: autoRegressionParameters.csv");

  //read header of file (unused)
  string csvLine;
  getline(f_autoRegressionParameters,csvLine);

  // read every line from the stream
  while(getline(f_autoRegressionParameters, csvLine)){

    std::istringstream csvStream(csvLine);

    string csvField1, csvField2, csvField3, csvField4, csvField5, csvField6, csvField7;

    // read every element from the line that is separated by commas
    getline(csvStream, csvField1, ',');
    getline(csvStream, csvField2, ',');
    getline(csvStream, csvField3, ',');
    getline(csvStream, csvField4, ',');
    getline(csvStream, csvField5, ',');
    getline(csvStream, csvField6, ',');
    getline(csvStream, csvField7, ',');

    istringstream csvNum1(csvField1), csvNum2(csvField2), csvNum3(csvField3), csvNum4(csvField4), csvNum5(csvField5), csvNum6(csvField6), csvNum7 (csvField7);
    
    int day;
    csvNum1 >> day;
    if (day >= _maximumDurationInDays)
      throw runtime_error ("EmpiricalInfection::init(): invalid day");
    csvNum2 >> _mu_beta1[day];
    csvNum3 >> _sigma_beta1[day];
    csvNum4 >> _mu_beta2[day];
    csvNum5 >> _sigma_beta2[day];
    csvNum6 >> _mu_beta3[day];
    csvNum7 >> _sigma_beta3[day];
  }  
  f_autoRegressionParameters.close();
}


// -----  non-static class members (Construction and destruction, checkpointing)  -----

/* Initialises a new infection by assigning the densities for the last 3 prepatent days
*/
EmpiricalInfection::EmpiricalInfection(double growthRateMultiplier) {
  //sample the parasite densities for the last 3 prepatent days
  //note that the lag decreases with time
  _laggedLogDensities[0]=sampleSubPatentValue(_alpha1,_mu1,log(_subPatentLimit));  
  _laggedLogDensities[1]=sampleSubPatentValue(_alpha2,_mu2,log(_subPatentLimit)); 
  _laggedLogDensities[2]=sampleSubPatentValue(_alpha3,_mu3,log(_subPatentLimit));  
  //only the immediately preceding value is modified by the growth rate multiplier
  _laggedLogDensities[0] += log(growthRateMultiplier); 
  _patentGrowthRateMultiplier = growthRateMultiplier;
  
  if (util::ModelOptions::option (util::INCLUDES_PK_PD))
      _proteome = PkPd::ProteomeInstance::newInfection();
}
EmpiricalInfection::~EmpiricalInfection() {
}


// -----  non-static class members (other functions)  -----

void EmpiricalInfection::setPatentGrowthRateMultiplier(double multiplier) {
  _patentGrowthRateMultiplier = multiplier;
}


// TODO: TS needs to improve documentation of this function.
const int EI_MAX_SAMPLES = 10;
bool EmpiricalInfection::updateDensity(int simulationTime, double survivalFactor) {
  //to avoid the formula for the linear predictor being excessively long we introduce L for the lagged densities
  # define L _laggedLogDensities
  
  int ageOfInfection = simulationTime - _startdate;	// age in days
  if (ageOfInfection > _maximumDurationInDays || !(L[0] > -999999.9))	// NOTE: second test is extremely unlikely to fail
    return true;	// cut-off point
  
  // constraints to ensure the density is defined and not exploding
  double upperLimitoflogDensity=log(_maximumPermittedAmplificationPerCycle*exp(L[1])/_inflationMean);
  double amplificationPerCycle;
  for (int tries0 = 0; tries0 < EI_MAX_SAMPLES; ++tries0) {
    double logDensity;
    for (int tries1 = 0; tries1 < EI_MAX_SAMPLES; ++tries1) {
      double b_1=gsl::rngGauss(_mu_beta1[ageOfInfection],_sigma_beta1[ageOfInfection]);
      double b_2=gsl::rngGauss(_mu_beta2[ageOfInfection],_sigma_beta2[ageOfInfection]);
      double b_3=gsl::rngGauss(_mu_beta3[ageOfInfection],_sigma_beta3[ageOfInfection]);
      double expectedlogDensity = b_1 * (L[0]+L[1]+L[2]) / 3
      + b_2 * (L[2]-L[0]) / 2
      + b_3 * (L[2]+L[0]-2*L[1]) / 4;
      
      //include sampling error
      logDensity=gsl::rngGauss(expectedlogDensity,sigma_noise(ageOfInfection));
      //include drug and immunity effects via growthRateMultiplier 
      logDensity += log(_patentGrowthRateMultiplier);
      
      if (logDensity <= upperLimitoflogDensity)	//got an acceptable density, we're done
	break;	// most of the time this should happen first try
    }
    if (!(logDensity <= upperLimitoflogDensity))	// in case all the above attempts fail, cap logDensity
      logDensity=upperLimitoflogDensity;
    
    _density= getInflatedDensity(logDensity); 
    
    _density *= survivalFactor;	// Apply drug and vaccine effects
    
    // Infections that get killed before they become patent:
    if ((ageOfInfection==0) && (_density < _subPatentLimit))
      _density=0.0;
    
    amplificationPerCycle=_density/exp(L[1]);
    if (_density >= 0.0 && amplificationPerCycle <= _maximumPermittedAmplificationPerCycle)
      break;	// We're done. Hopefully usually with the first try.
  }
  if (!(_density >= 0.0 && amplificationPerCycle <= _maximumPermittedAmplificationPerCycle))	// in case the above tries fail
    _density = _maximumPermittedAmplificationPerCycle*exp(L[1]);
  
  _laggedLogDensities[2]=_laggedLogDensities[1];
  _laggedLogDensities[1]=_laggedLogDensities[0];
  _laggedLogDensities[0]=log(_density);
  
  _cumulativeExposureJ += Global::interval * _density;
  
  // Note: here use a positive test for survival, since if _density became an NaN tests against it will return false:
  if (_density*_overallMultiplier > _extinctionLevel)
    return false;	// Still parasites; infection didn't go extinct
  else
    return true;	// parasites are extinct; infection will be removed from model
# undef L
}

double EmpiricalInfection::sampleSubPatentValue(double alpha, double mu, double upperBound){
  double beta=alpha*(1-mu)/mu;
  double nonInflatedValue=upperBound+log(gsl::rngBeta(alpha, beta));
  double inflatedValue;
  int tries=0;
  do {
    inflatedValue=getInflatedDensity(nonInflatedValue);
    tries++;
  } while ((inflatedValue>upperBound) && tries<10);
  if (inflatedValue>upperBound) inflatedValue=upperBound;
  return inflatedValue;
}

double EmpiricalInfection::samplePatentValue(double mu, double sigma, double lowerBound){
  double returnValue;
  do {
    double nonInflatedValue=gsl::rngGauss(mu, sigma);
    returnValue=getInflatedDensity(nonInflatedValue);
  } while (returnValue<lowerBound);
  return returnValue;
}

double EmpiricalInfection::sigma_noise(int ageOfInfection) {
  return _sigma0_res+_sigmat_res*((double)ageOfInfection);
}

double EmpiricalInfection::getInflatedDensity(double nonInflatedDensity){
  double inflatedLogDensity=log(_inflationMean)+gsl::rngGauss(nonInflatedDensity,sqrt(_inflationVariance));
  return exp(inflatedLogDensity);
}

void EmpiricalInfection::overrideInflationFactors(double inflationMean, double inflationVariance, double extinctionLevel, double overallMultiplier){
  _inflationVariance=inflationVariance;
  _inflationMean=inflationMean;
  _extinctionLevel=extinctionLevel;
  _overallMultiplier=overallMultiplier;
  _subPatentLimit=10.0/_overallMultiplier;
}


// -----  checkpointing  -----

void EmpiricalInfection::checkpoint (istream& stream) {
    Infection::checkpoint (stream);
    _laggedLogDensities[0] & stream;
    _laggedLogDensities[1] & stream;
    _laggedLogDensities[2] & stream;
    _patentGrowthRateMultiplier & stream;
}
void EmpiricalInfection::checkpoint (ostream& stream) {
    Infection::checkpoint (stream);
    _laggedLogDensities[0] & stream;
    _laggedLogDensities[1] & stream;
    _laggedLogDensities[2] & stream;
    _patentGrowthRateMultiplier & stream;
}

} }