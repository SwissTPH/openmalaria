/*
 This file is part of OpenMalaria.

 Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine

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

#include "WithinHost/DescriptiveInfection.h"
#include "Host/intervention.h"
#include "inputData.h"
#include "util/random.h"
#include "util/CommandLine.hpp"
#include "util/ModelOptions.hpp"
#include "util/errors.hpp"
#include <algorithm>
#include <sstream>
#include <string.h>
#include <stdexcept>

namespace OM { namespace WithinHost {
    using namespace util;
    //static (class) variables

int DescriptiveInfection::latentp;
double DescriptiveInfection::meanLogParasiteCount[maxDur*maxDur];
double DescriptiveInfection::sigma0sq;
double DescriptiveInfection::xNuStar;


// -----  static init/clear -----

void DescriptiveInfection::initParameters (){
  if (Global::interval != 5)
    throw util::xml_scenario_error ("DescriptiveInfection only supports using an interval of 5");
  if (util::ModelOptions::option (util::INCLUDES_PK_PD))
      throw util::xml_scenario_error ("INCLUDES_PK_PD is incompatible with the old within-host model");
  
  latentp=InputData().getModel().getParameters().getLatentp();
  sigma0sq=InputData.getParameter(Params::SIGMA0_SQ);
  xNuStar=InputData.getParameter(Params::X_NU_STAR);
  //File name of file with empirical parasite densities.
  string densities_filename;
  densities_filename = util::CommandLine::lookupResource ("densities.csv");

  fstream f_MTherapyDensities(densities_filename.c_str(),ios::in);

  if (!f_MTherapyDensities.is_open()){
    //If File cannot be accessed:
    cerr << "file not found: densities.csv" << endl;
    exit(-1);
  }

  //read header of file (unused)
  string csvLine;
  getline(f_MTherapyDensities,csvLine);

  // read every line from the stream
  while(getline(f_MTherapyDensities, csvLine)){
    //Empirical description of single Malaria infections in naive individuals
    //counter variables, i stands for 5 day time interval, j for duration of infection
    int i;
    int j;
    double meanlogdens;

    std::istringstream csvStream(csvLine);

    string csvField1, csvField2, csvField3;

    // read every element from the line that is seperated by commas
    getline(csvStream, csvField1, ',');
    getline(csvStream, csvField2, ',');
    getline(csvStream, csvField3, ',');

    istringstream csvNum1(csvField1), csvNum2(csvField2), csvNum3(csvField3);
   
    csvNum1 >> i;
    csvNum2 >> j;
    csvNum3 >> meanlogdens;

    //fill initial matrix
    meanLogParasiteCount[i-1+(j-1)*maxDur]=meanlogdens;
    //fill also the triangle that will not be used (to ensure everything is initialised)
    if (j!=i){
      meanLogParasiteCount[j-1+(i-1)*maxDur]=0.0;
    }

  }

  f_MTherapyDensities.close();
}

void DescriptiveInfection::clearParameters () {}


// -----  non-static init/destruction  -----

DescriptiveInfection::DescriptiveInfection () :
    Infection(0xFFFFFFFF),
    _startdate(Global::simulationTime)
{
    _duration=infectionDuration();
}

DescriptiveInfection::~DescriptiveInfection() {
}


// -----  other  -----

int DescriptiveInfection::infectionDuration(){
    double meanlogdur=5.1300001144409179688;
    //Std of the logduration
    double sdlogdur=0.80000001192092895508;
    double dur=random::log_normal(meanlogdur, sdlogdur);
    return (1+(int)floor(dur)) / Global::interval;
}

void DescriptiveInfection::determineDensities(double ageInYears, double cumulativeh, double cumulativeY, double &timeStepMaxDensity, double innateImmSurvFact, double BSVEfficacy)
{
  //Age of infection. (Blood stage infection starts latentp intervals later than inoculation ?)
  int infage=1+Global::simulationTime-_startdate-latentp;
  if ( infage >  0) {
    if ( infage <=  maxDur) {
      int iduration=_duration;
      if ( iduration >  maxDur)
	iduration=maxDur;
      
      _density=exp(meanLogParasiteCount[infage - 1 + (iduration - 1)*maxDur]);
    } else {
      _density=exp(meanLogParasiteCount[maxDur - 1 + (maxDur - 1)*maxDur]);
    }
    if (_density < 1.0)
      _density=1.0;
    
    /*
    The expected parasite density in the non naive host. 
    As regards the second term in AJTM p.9 eq. 9, in published and current implementations Dx is zero.
    */
    _density = exp(log(_density) * immunitySurvivalFactor(ageInYears, cumulativeh, cumulativeY));
    //Perturb _density using a lognormal 
    double varlog = sigma0sq / (1.0 + (cumulativeh / xNuStar));
    double stdlog = sqrt(varlog);
    /*
    This code samples from a log normal distribution with mean equal to the predicted density
    n.b. AJTM p.9 eq 9 implies that we sample the log of the density from a normal with mean equal to
    the log of the predicted density.  If we really did the latter then this bias correction is not needed.
    */
    double meanlog = log(_density) - stdlog*stdlog / 2.0;
    timeStepMaxDensity = 0.0;
    if (stdlog > 0.0000001) {
      if (Global::interval > 1) {
	/*
	sample the maximum density over the T-1 remaining days in the
	time interval, (where T is the duration of the time interval)
	*/
	double normp = pow(random::uniform_01(), 1.0 / (Global::interval-1));
	/*
	To mimic sampling T-1 repeated values, we transform the sampling
	distribution and use only one sampled value, which has the sampling
	distribution of the maximum of T-1 values sampled from a uniform.
	The probability density function of this sampled random var is distributed
	according to a skewed distribution (defined in [0,1]) where the
	exponent (1/(T-1)) arises because each of T-1 sampled
	values would have this probability of being the maximum. 
	*/
	timeStepMaxDensity = random::sampleFromLogNormal(normp, meanlog, stdlog);
      }
      //calculate the expected density on the day of sampling
      _density = random::sampleFromLogNormal(random::uniform_01(), meanlog, stdlog);
      timeStepMaxDensity = std::max(_density, timeStepMaxDensity);
    }
    if (_density > maxDens || timeStepMaxDensity > maxDens) {
      cerr << "MD lim" << endl;
      _density = maxDens;
      timeStepMaxDensity = _density;
    }
  }
  else {
    _density = 0.0;
  }
  
  //Compute the proportion of parasites remaining after innate blood stage effect
  _density *= innateImmSurvFact;
  /* MAX_DENS_BUG: Possibly a better model version ensuring that the effect of
   * variation in innate immunity is reflected in case incidence would have the
   * following: */
  if (util::ModelOptions::option (util::INNATE_MAX_DENS))
      timeStepMaxDensity *= innateImmSurvFact;
  
  //Include here the effect of blood stage vaccination
  if (Host::Vaccine::BSV.active) {
    double factor = 1.0 - BSVEfficacy;
    _density *= factor;
    timeStepMaxDensity *= factor;
  }
}

//Note: would make sense is this was also part of determineDensities, but can't really be without changing order of other logic.
void DescriptiveInfection::determineDensityFinal () {
  _density = std::min(maxDens, _density);
  _cumulativeExposureJ += Global::interval * _density;
}


DescriptiveInfection::DescriptiveInfection (istream& stream) :
    Infection(stream)
{
    _startdate & stream;
    _duration & stream;
}
void DescriptiveInfection::checkpoint (ostream& stream) {
    Infection::checkpoint (stream);
    _startdate & stream;
    _duration & stream;
}

} }