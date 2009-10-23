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
#include "inputData.h"
#include "util/gsl.h"
#include <algorithm>
#include <sstream>
#include <string.h>

//static (class) variables

int DescriptiveInfection::latentp;
double DescriptiveInfection::meanLogParasiteCount[maxDur*maxDur];
double DescriptiveInfection::alpha_m;
double DescriptiveInfection::decayM;
double DescriptiveInfection::sigma0sq;
double DescriptiveInfection::xNuStar;


// -----  static init/clear -----

void DescriptiveInfection::initParameters (){
  if (Global::interval != 5)
    throw domain_error ("DescriptiveInfection only supports using an interval of 5");
  
  latentp=get_latentp();
  cumulativeYstar=(float)getParameter(Params::CUMULATIVE_Y_STAR);
  cumulativeHstar=(float)getParameter(Params::CUMULATIVE_H_STAR);
  alpha_m=1-exp(-getParameter(Params::NEG_LOG_ONE_MINUS_ALPHA_M));
  decayM=getParameter(Params::DECAY_M);
  sigma0sq=getParameter(Params::SIGMA0_SQ);
  xNuStar=getParameter(Params::X_NU_STAR);
  //File name of file with empirical parasite densities.
  string densities_filename;
  densities_filename = Global::lookupResource ("densities.csv");

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

DescriptiveInfection::DescriptiveInfection(int simulationTime) :
  Infection(simulationTime)
{
    //Initialize current infection data
    _duration=infectionDuration();
    _cumulativeExposureJ=0.0;
    
    if (Global::modelVersion & INCLUDES_PK_PD)
      _proteome = ProteomeInstance::newInfection();
    else
      _proteome = NULL;
}

DescriptiveInfection::~DescriptiveInfection() {
}

DescriptiveInfection::DescriptiveInfection (istream& in) :
  Infection (in)
{
  in >> _duration;
  in >> _cumulativeExposureJ; 
}

void DescriptiveInfection::write (ostream& out) const {
  Infection::write (out);
  out << _duration << endl; 
  out << _cumulativeExposureJ << endl; 
}


// -----  other  -----

int DescriptiveInfection::getEndDate(){
  return _startdate+_duration/Global::interval;
}

int DescriptiveInfection::infectionDuration(){
    double meanlogdur=5.1300001144409179688;
    //Std of the logduration
    double sdlogdur=0.80000001192092895508;
    double dur=gsl::rngLogNormal(meanlogdur, sdlogdur);
    return 1+(int)floor(dur);
}

void DescriptiveInfection::determineDensities(int simulationTime, double cumulativeY, double ageYears, double cumulativeh, double &timeStepMaxDensity)
{
  //Age of infection. (Blood stage infection starts latentp intervals later than inoculation ?)
  int infage=1+simulationTime-_startdate-latentp;
  if ( infage >  0) {
    if ( infage <=  maxDur) {
      int iduration=_duration/Global::interval;
      if ( iduration >  maxDur)
	iduration=maxDur;
      
      _density=exp(meanLogParasiteCount[infage - 1 + (iduration - 1)*maxDur]);
    } else {
      _density=exp(meanLogParasiteCount[maxDur - 1 + (maxDur - 1)*maxDur]);
    }
    if (_density < 1.0)
      _density=1.0;
    
    //effect of cumulative Parasite density (named Dy in AJTM)
    double dY;
    //effect of number of infections experienced since birth (named Dh in AJTM)
    double dH;
    //effect of age-dependent maternal immunity (named Dm in AJTM)
    double dA;
    
    if (cumulativeh <= 1.0) {
      dY=1.0;
      dH=1.0;
    } else {
      dH=1.0 / (1.0 + (cumulativeh-1.0) / cumulativeHstar);
      //TODO: compare this with the asex paper
      dY=1.0 / (1.0 + (cumulativeY-_cumulativeExposureJ) / cumulativeYstar);
    }
    dA = 1.0 - alpha_m * exp(-decayM * ageYears);
    double survival = std::min(dY*dH*dA, 1.0);
    /*
    The expected parasite density in the non naive host. 
    As regards the second term in AJTM p.9 eq. 9, in published and current implementations Dx is zero.
    */
    _density = exp(log(_density) * survival);
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
	double normp = pow(gsl::rngUniform(), 1.0 / (Global::interval-1));
	/*
	To mimic sampling T-1 repeated values, we transform the sampling
	distribution and use only one sampled value, which has the sampling
	distribution of the maximum of T-1 values sampled from a uniform.
	The probability density function of this sampled random var is distributed
	according to a skewed distribution (defined in [0,1]) where the
	exponent (1/(T-1)) arises because each of T-1 sampled
	values would have this probability of being the maximum. 
	*/
	timeStepMaxDensity = gsl::sampleFromLogNormal(normp, meanlog, stdlog);
      }
      //calculate the expected density on the day of sampling
      _density = gsl::sampleFromLogNormal(gsl::rngUniform(), meanlog, stdlog);
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
}
