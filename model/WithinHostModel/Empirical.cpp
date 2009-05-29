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

#include "WithinHostModel/Empirical.h"
#include "GSLWrapper.h"
#include <sstream>
#include <fstream>


EmpiricalWithinHostModel::EmpiricalWithinHostModel(){
  // alpha1 corresponds to 1 day before first patent, alpha2 2 days before first patent etc.
  _lambda=-1.1;
  _alpha1=	-2.909;
	_alpha2=	-0.9317;	
  _alpha3=	2.72;	
  _sigma_alpha1=0.001091;
  _sigma_alpha2=0.005054;
  _sigma_alpha3=2.634;
  _sigma0_res=0.02334;
  _sigmat_res=0.002157;
  _maximumPermittedAmplificationPerCycle=1000.0;
  double _subPatentLimit=10.0;
  const int MAX_LENGTH = 1000;
  char autoRegressionParameters [MAX_LENGTH];
  strcpy(autoRegressionParameters,"autoRegressionParameters.csv");
  fstream f_autoRegressionParameters(autoRegressionParameters,ios::in);
  int day;
  if (!f_autoRegressionParameters.is_open()){
    //If File cannot be accessed:
    cerr << "file not found: autoRegressionParameters.csv" << endl;
    exit(-1);
  }

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
   
    csvNum1 >> day;
    csvNum2 >> _mu_beta1[day];
    csvNum3 >> _sigma_beta1[day];
    csvNum4 >> _mu_beta2[day];
    csvNum5 >> _sigma_beta2[day];
    csvNum6 >> _mu_beta3[day];
    csvNum7 >> _sigma_beta3[day];
  }  
  f_autoRegressionParameters.close();
  
}

void EmpiricalWithinHostModel::initialiseInfection(double * transformedLaggedDensities){
//sample the parasite densities for the last 3 prepatent days
//note that the lag decreases with time
transformedLaggedDensities[0]=sampleSubPatentValue(_alpha1,_sigma_alpha1,boxCoxTransform(_subPatentLimit));  
transformedLaggedDensities[1]=sampleSubPatentValue(_alpha2,_sigma_alpha2,boxCoxTransform(_subPatentLimit)); 
transformedLaggedDensities[2]=sampleSubPatentValue(_alpha3,_sigma_alpha3,boxCoxTransform(_subPatentLimit));  
return;  
}

double EmpiricalWithinHostModel::getNewDensity(double * transformedLaggedDensities, int ageOfInfection){
double newDensity=-9.99;  
double transformedInflatedDensity=-9999999.99;
if (transformedLaggedDensities[0]>-999999.9) {
//to avoid the formula for the linear predictor being excessively long we introduce L for the lagged densities
  double L[3];
  for (int i=0;i<3;i++) L[i]=transformedLaggedDensities[i];
// constraints to ensure the density is defined and not exploding
  double upperLimitofTransformedDensity=boxCoxTransform(_maximumPermittedAmplificationPerCycle*inverseBoxCoxTransform(L[1])/_inflationMean);
  double amplificationPerCycle=999999.9;
  int tries0=0;
  while (((newDensity <0) || (amplificationPerCycle > _maximumPermittedAmplificationPerCycle)) && (tries0<10)){
    int tries1=0;
    double transformedDensity=9999.9;
    while ((transformedDensity>upperLimitofTransformedDensity) && (tries1<10)) {
      double b_1=W_GAUSS(_mu_beta1[ageOfInfection],_sigma_beta1[ageOfInfection]);
      double b_2=W_GAUSS(_mu_beta2[ageOfInfection],_sigma_beta2[ageOfInfection]);
      double b_3=W_GAUSS(_mu_beta3[ageOfInfection],_sigma_beta3[ageOfInfection]);
      double expectedTransformedDensity=b_1*(L[0]+L[1]+L[2])/3+b_2*(L[2]-L[0])/2+b_3*(L[2]+L[0]-2*L[1])/4;
//include sampling error
      transformedDensity=W_GAUSS(expectedTransformedDensity,sigma_noise(ageOfInfection));
      tries1++;
    }
    if (tries1 > 9) transformedDensity=upperLimitofTransformedDensity;
    newDensity= getInflatedDensity(transformedDensity); 
    if ((ageOfInfection==0) && (newDensity < _subPatentLimit)) newDensity=-9.9; 
    tries0++;
    if (tries0 > 9) newDensity=_maximumPermittedAmplificationPerCycle*inverseBoxCoxTransform(L[1]);
    transformedInflatedDensity=boxCoxTransform(newDensity);
    amplificationPerCycle=newDensity/inverseBoxCoxTransform(L[1]);
  }
}
transformedLaggedDensities[2]=transformedLaggedDensities[1];
transformedLaggedDensities[1]=transformedLaggedDensities[0];
transformedLaggedDensities[0]=transformedInflatedDensity;
if (newDensity*_overallMultiplier< _extinctionLevel) {
  transformedLaggedDensities[0]=-9999999.99;
  newDensity=-9.99;
}
return newDensity*_overallMultiplier;
}

double EmpiricalWithinHostModel::sampleSubPatentValue(double mu, double sigma, double upperBound){
  double nonInflatedValue;
  int tries=0;
  do {
    nonInflatedValue=W_GAUSS(mu, sigma);
    tries++;
  } while ((nonInflatedValue>upperBound) && tries<10);
  if (nonInflatedValue>upperBound) nonInflatedValue=upperBound;
  double inflatedValue;
  tries=0;
  do {
    inflatedValue=getInflatedDensity(nonInflatedValue);
    tries++;
  } while ((inflatedValue>upperBound) && tries<10);
  if (inflatedValue>upperBound) inflatedValue=upperBound;
  return inflatedValue;
}

double EmpiricalWithinHostModel::samplePatentValue(double mu, double sigma, double lowerBound){
  double returnValue;
  do {
    double nonInflatedValue=W_GAUSS(mu, sigma);
    returnValue=getInflatedDensity(nonInflatedValue);
  } while (returnValue<lowerBound);
  return returnValue;
}

double EmpiricalWithinHostModel::inverseBoxCoxTransform(double transformedValue){
return exp(log(_lambda*transformedValue+1.0)/_lambda);
}

double EmpiricalWithinHostModel::boxCoxTransform(double untransformedValue){
return (exp(_lambda*log(untransformedValue))-1.0)/_lambda;
}

void EmpiricalWithinHostModel::setInflationFactors(double inflationMean, double inflationVariance, double extinctionLevel, double overallMultiplier){
  _inflationVariance=inflationVariance;
  _inflationMean=inflationMean;
  _extinctionLevel=extinctionLevel;
  _overallMultiplier=overallMultiplier;
  _subPatentLimit=10.0/overallMultiplier;
}

double EmpiricalWithinHostModel::sigma_noise(int ageOfInfection) {
  return _sigma0_res+_sigmat_res*((double)ageOfInfection);
}

double EmpiricalWithinHostModel::getInflatedDensity(double nonInflatedDensity){
  double logBackTransformedDensity= log(inverseBoxCoxTransform(nonInflatedDensity));  
  double inflatedLogDensity=log(_inflationMean)+W_GAUSS(logBackTransformedDensity,sqrt(_inflationVariance));
return exp(inflatedLogDensity);
}