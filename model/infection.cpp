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

#include "infection.h"
#include "inputData.h"
#include "GSLWrapper.h"
#include <iostream>
#include <algorithm>
#include <sstream>
#include <string.h>

#include "boinc_bridge.h"

//static (class) variables

float Infection::cumulativeYstar;
float Infection::cumulativeHstar;
double Infection::meanLogParasiteCount[maxDur*maxDur];
double Infection::alpha_m;
double Infection::decayM;
double Infection::sigma0sq;
double Infection::xNuStar;

Infection::Infection() {
}

Infection::~Infection() {
  //if ( isOptionIncluded(modelVersion, includesPKPD)) {
  //  delete _proteome;
  //}
}

void Infection::initParameters(){
  //Empirical description of single Malaria infections in naive individuals
  //counter variables, i stands for 5 day time interval, j for duration of infection
  int i;
  int j;
  double meanlogdens;
  cumulativeYstar=(float)get_parameter(7);
  cumulativeHstar=(float)get_parameter(8);
  alpha_m=1-exp(-get_parameter(9));
  decayM=get_parameter(10);
  sigma0sq=get_parameter(11);
  xNuStar=get_parameter(12);
  //File name of file with empirical parasite densities.
  string densities_filename;
  int retval = boinc_resolve_filename_s("densities.csv",densities_filename);
  if (retval){
    std::cerr << "APP. boinc_resolve_filename failed \n";
    boinc_finish(retval);
  }

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

Infection::Infection(int lastSPdose, int simulationTime){
    int genotypeCounter;
    double uniformRandomVariable;
    double lowerIntervalBound;
    double upperIntervalBound;
    //Initialize current infection data
    _startdate=simulationTime;
    _density=0.0;
    _duration=infectionDuration();
    _cumulativeExposureJ=0.0;
    _gType.ID=0;
    if ( IPT) {
        uniformRandomVariable=(W_UNIFORM());
        lowerIntervalBound=0.0;
        upperIntervalBound=genotypeFreq[1 - 1];
        //This Loop assigns the infection a genotype according to its frequency
        for ( genotypeCounter=1;genotypeCounter<=numberOfGenoTypes; genotypeCounter++){
            if ( uniformRandomVariable > lowerIntervalBound && uniformRandomVariable < upperIntervalBound){
                _gType.ID=genotypeCounter;
            }
            lowerIntervalBound=upperIntervalBound;
            /*
            The loop should never come into this else-part (exits before), so the if statement is not necessary.
            For safety reason we do it nevertheless
            */
            if ( genotypeCounter !=  numberOfGenoTypes) {
                upperIntervalBound=upperIntervalBound+genotypeFreq[genotypeCounter+1 - 1];
            }
            else {
              upperIntervalBound=1.0;
              //write(0,*) 'should not be here'
            }
        }
        /*
        The attenuation effect of SP is only effective during a certain time-window for certain IPTi models
        If t(=now) lies within this time window, SPattenuate is true, false otherwise.
        The time window starts after the prophylactic period ended (during the prophylactic
        period infections are cleared) and ends genotypeTolPeriod(iTemp%iData%gType%ID) time steps later.
        */
        if (simulationTime-lastSPdose > genotypeProph[_gType.ID-1] && simulationTime-lastSPdose <= genotypeProph[_gType.ID-1]+genotypeTolPeriod[_gType.ID-1]){
          _SPattenuate=true;
        }
    }
    //This should probably be inside an IF
    if ( isOptionIncluded(modelVersion, includesPKPD)) {
      _proteome = ProteomeManager::getManager()->getInfection();
    }
    else {
      _proteome = 0;
    }
}

void Infection::writeInfectionToFile(fstream& funit){
  funit << *this;
}

int Infection::getEndDate(){
return _startdate+_duration/interval;
}

double Infection::determineWithinHostDensity(){
    // TODO: this is the insertion point for the within host model

    double density;
    double valdetermineWithinHostDensity;
    density=std::max(0.025, _density);
    /*
     if the density gets to be < 1 parasite per host then clear infections
     is called by making the duration negative. TODO: the current test version
     does not give a reasonable model for densities so we do not allow
     the within-host model to set the duration.  Instead we allow the original model to determine duration
     An alternative would be to set initial duration to be large positive to override the original model
    */
    if ( density <  0.02) {
        _duration=-99;
    }
    valdetermineWithinHostDensity=mymodf(density*8.0, 20000.0);
    return valdetermineWithinHostDensity;
}

double sampleFromLogNormal(double normp, double meanlog, double stdlog){
    /*
    
    Used for performance reasons. Calling GLS LOG_NORMAL 5 times is 50% slower.
    
    */

    double zval;
    double valsampleFromLogNormal;
    zval=W_UGAUSS_PINV(normp);
    /*
    Why not zval=W_UGAUSS?
    where normp is distributed uniformly over [0,1],
    zval is distributed like a standard normal distribution
    where normp has been transformed by raising to the power of 1/(T-1) 
    zval is distributed like a uniform gauss	times 4* F(x,0,1)^3, where F(x,0,1) ist the cummulative
    distr. function of a uniform gauss
    */
    valsampleFromLogNormal=exp(meanlog+stdlog*((float)zval));
    return valsampleFromLogNormal;
}

int Infection::infectionDuration(){
    double meanlogdur=5.1300001144409179688;
    //Std of the logduration
    double sdlogdur=0.80000001192092895508;
    int valinfectionDuration;
    double dur=W_LOGNORMAL(meanlogdur, sdlogdur);
    valinfectionDuration=1+(int)floor(dur);
    return valinfectionDuration;
}

ProteomeInstance* Infection::getProteome() const {
  return _proteome;
}

ostream& operator<<(ostream& out, const Infection& infection){
  
  out << infection._duration << endl; 
  out << infection._startdate << endl; 
  out << infection._density << endl; 
  out << infection._cumulativeExposureJ << endl; 
  out << infection._gType.ID << endl; 
  if ( isOptionIncluded(modelVersion, includesPKPD)) {
    out << infection._proteome->getProteomeID() << endl; 
  }
  return out << boolalpha << infection._SPattenuate << endl; 

}

istream& operator>>(istream& in, Infection& infection){
  int proteomeID;
  in >> infection._duration; 
  in >> infection._startdate; 
  in >> infection._density; 
  in >> infection._cumulativeExposureJ; 
  in >> infection._gType.ID; 
  if ( isOptionIncluded(modelVersion, includesPKPD)) {
    in >> proteomeID; 
    infection._proteome = ProteomeManager::getManager()->getProteome(proteomeID);
  }
  in >> boolalpha >> infection._SPattenuate; 

  return in;

}

void Infection::determineDensities(int simulationTime, double cumulativeY, double ageyears, double cumulativeh, double *timeStepMaxDensity){

    int iduration;
    int infage;
    double survival;
    double normp;
    double y;
    double logy;
    double stdlog;
    double meanlog;
    double varlog;
    //effect of cumulative Parasite density (named Dy in AJTM)
    double dY;
    //effect of number of infections experienced since birth (named Dh in AJTM)
    double dH;
    //effect of age-dependent maternal immunity (named Dm in AJTM)
    double dA;
    //Age of infection. (Blood stage infection starts latentp intervals later than inoculation ?)
    infage=1+simulationTime-_startdate-nearbyint(latentp);
    if ( infage >  0) {
        iduration=_duration/interval;
        if ( iduration >  maxDur) {
            iduration=maxDur;
        }
        if ( infage <=  maxDur) {
            y=(float)exp(meanLogParasiteCount[infage - 1 + (iduration - 1)*maxDur]);
        }
        else {
            y=(float)exp(meanLogParasiteCount[maxDur - 1 + (maxDur - 1)*maxDur]);
        }
        if ( y <  1.0) {
            y=1.0;
        }
        if ( cumulativeh <=  1.0) {
            dY=1;
            dH=1;
        }
        else {
            dH=1/(1+(cumulativeh-1.0)/cumulativeHstar);
            //TODO: compare this with the asex paper
            dY=1/(1+(cumulativeY-_cumulativeExposureJ)/cumulativeYstar);
        }
        //Can this happen or is it just for security ?
        if ( ageyears <=  0.0) {
            dA=1-alpha_m;
        }
        else {
            dA=1-alpha_m*exp(-decayM*ageyears);
        }
        survival=dY*dH*dA;
        survival=std::min(survival, 1.0);
        logy=log(y)*(survival);
        /*
        The expected parasite density in the non naive host. 
        As regards the second term in AJTM p.9 eq. 9, in published and current implementations Dx is zero.
        */
        y=exp(logy);
        //Perturb y using a lognormal 
        varlog=sigma0sq/(1+(cumulativeh/xNuStar));
        stdlog=sqrt(varlog);
        /*
        This code samples from a log normal distribution with mean equal to the predicted density
        n.b. AJTM p.9 eq 9 implies that we sample the log of the density from a normal with mean equal to
        the log of the predicted density.  If we really did the latter then this bias correction is not needed.
        */
        meanlog=log(y)-stdlog*stdlog/2.0;
        *timeStepMaxDensity =0.0;
        if ( stdlog >  0.0000001) {
            if ( interval >  1) {
                normp=W_UNIFORM();
                /*
                sample the maximum density over the T-1 remaining days in the
                time interval, (where T is the duration of the time interval)
                */
                normp=pow(normp, 1.0*1/(interval-1));
                /*
                To mimic sampling T-1 repeated values, we transform the sampling
                distribution and use only one sampled value, which has the sampling
                distribution of the maximum of T-1 values sampled from a uniform.
                The probability density function of this sampled random var is distributed
                according to a skewed distribution (defined in [0,1]) where the
                exponent (1/(T-1)) arises because each of T-1 sampled
                values would have this probability of being the maximum. 
                */
                *timeStepMaxDensity =sampleFromLogNormal(normp, meanlog, stdlog);
            }
            //calculate the expected density on the day of sampling
            normp=W_UNIFORM();
            y=(float)sampleFromLogNormal(normp, meanlog, stdlog);
            *timeStepMaxDensity = std::max((double) y, *timeStepMaxDensity);
        }
        if (( y >  maxDens) || ( *timeStepMaxDensity >  (double)maxDens)) {
          cout << "MD lim" << endl;
          y=maxDens;
          *timeStepMaxDensity =(double) y;
        }
        _density = y;
    }
    else {
        _density =0.0;
    }
}
