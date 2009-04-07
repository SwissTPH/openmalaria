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

#include "DescriptiveInfection.h"
#include "inputData.h"
#include "GSLWrapper.h"
#include <iostream>
#include <algorithm>
#include <sstream>
#include <string.h>

#include "boincWrapper.h"

//static (class) variables

double DescriptiveInfection::meanLogParasiteCount[maxDur*maxDur];
double DescriptiveInfection::alpha_m;
double DescriptiveInfection::decayM;
double DescriptiveInfection::sigma0sq;
double DescriptiveInfection::xNuStar;

// static IPT variables (possibly move to a subclass)
bool DescriptiveInfection::IPT;
int DescriptiveInfection::numberOfGenoTypes;
double *DescriptiveInfection::genotypeFreq;
int *DescriptiveInfection::genotypeProph;
int *DescriptiveInfection::genotypeTolPeriod;
double *DescriptiveInfection::genotypeACR;
double *DescriptiveInfection::genotypeAtten;


// -----  static init/clear -----

void DescriptiveInfection::initParameters (){
  cumulativeYstar=(float)getParameter(Params::CUMULATIVE_Y_STAR);
  cumulativeHstar=(float)getParameter(Params::CUMULATIVE_H_STAR);
  alpha_m=1-exp(-getParameter(Params::NEG_LOG_ONE_MINUS_ALPHA_M));
  decayM=getParameter(Params::DECAY_M);
  sigma0sq=getParameter(Params::SIGMA0_SQ);
  xNuStar=getParameter(Params::X_NU_STAR);
  //File name of file with empirical parasite densities.
  string densities_filename;
  densities_filename = BoincWrapper::resolveFile("densities.csv");

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
  
  
  // IPT-specific intitialisation
  
  const Interventions& xmlInterventions = getInterventions();
  IPT = xmlInterventions.getIptiDescription().present();
  if (!IPT)
    return;
  
  // --- IptiDescription begin ---
  const IptDescription& xmlIPTI = xmlInterventions.getIptiDescription().get();
  
  const IptDescription::InfGenotypeSequence& genotypes = xmlIPTI.getInfGenotype();
  numberOfGenoTypes = genotypes.size();
  genotypeFreq	= (double*)malloc(((numberOfGenoTypes))*sizeof(double));
  genotypeACR	= (double*)malloc(((numberOfGenoTypes))*sizeof(double));
  genotypeProph	= (int*)malloc(((numberOfGenoTypes))*sizeof(int));
  genotypeTolPeriod = (int*)malloc(((numberOfGenoTypes))*sizeof(int));
  genotypeAtten	= (double*)malloc(((numberOfGenoTypes))*sizeof(double));
  
  size_t i = 0;
  for (IptDescription::InfGenotypeConstIterator it = genotypes.begin(); it != genotypes.end(); ++it, ++i) {
    genotypeFreq[i]	= it->getFreq();
    genotypeACR[i]	= it->getACR();
    genotypeProph[i]	= it->getProph();
    genotypeTolPeriod[i]= it->getTolPeriod();
    genotypeAtten[i]	= it->getAtten();
  }
  // --- IptiDescription end ---
}

void DescriptiveInfection::clearParameters () {
  if (!IPT)
    return;
  
  free(genotypeFreq);
  free(genotypeACR);
  free(genotypeProph);
  free(genotypeTolPeriod);
  free(genotypeAtten);
}


// -----  non-static init/destruction  -----

DescriptiveInfection::DescriptiveInfection(int lastSPdose, int simulationTime){
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
    if (IPT) {
        uniformRandomVariable=(W_UNIFORM());
        lowerIntervalBound=0.0;
        upperIntervalBound=genotypeFreq[0];
        //This Loop assigns the infection a genotype according to its frequency
        for (genotypeCounter=1; genotypeCounter<=numberOfGenoTypes; genotypeCounter++){
            if (uniformRandomVariable > lowerIntervalBound &&
                uniformRandomVariable < upperIntervalBound){
                _gType.ID=genotypeCounter;
            }
            lowerIntervalBound=upperIntervalBound;
            /*
            The loop should never come into this else-part (exits before), so the if statement is not necessary.
            For safety reason we do it nevertheless
            */
            if ( genotypeCounter !=  numberOfGenoTypes) {
              upperIntervalBound = upperIntervalBound + genotypeFreq[genotypeCounter];
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
        if (simulationTime-lastSPdose > genotypeProph[_gType.ID-1] &&
            simulationTime-lastSPdose <= genotypeProph[_gType.ID-1] + genotypeTolPeriod[_gType.ID-1]){
          _SPattenuate=true;
        }
    }
    //This should probably be inside an IF
    if (Global::modelVersion & INCLUDES_PK_PD) {
      _proteome = ProteomeManager::getManager()->getInfection();
    }
    else {
      _proteome = 0;
    }
}

DescriptiveInfection::~DescriptiveInfection() {
}
void DescriptiveInfection::destroy() {
  //if (modelVersion & INCLUDES_PK_PD) {
  //  delete _proteome;
  //}
}


// -----  other  -----

void DescriptiveInfection::writeInfectionToFile(fstream& funit){
  write(funit);
}

int DescriptiveInfection::getEndDate(){
  return _startdate+_duration/Global::interval;
}

int DescriptiveInfection::infectionDuration(){
    double meanlogdur=5.1300001144409179688;
    //Std of the logduration
    double sdlogdur=0.80000001192092895508;
    int valinfectionDuration;
    double dur=W_LOGNORMAL(meanlogdur, sdlogdur);
    valinfectionDuration=1+(int)floor(dur);
    return valinfectionDuration;
}

void DescriptiveInfection::write (ostream& out) const {
  out << _duration << endl; 
  out << _startdate << endl; 
  out << _density << endl; 
  out << _cumulativeExposureJ << endl; 
  out << _gType.ID << endl; 
  if (Global::modelVersion & INCLUDES_PK_PD) {
    out << _proteome->getProteomeID() << endl; 
  }
  out << boolalpha << _SPattenuate << endl; 
}

DescriptiveInfection::DescriptiveInfection (istream& in) {
  int proteomeID;
  in >> _duration; 
  in >> _startdate; 
  in >> _density; 
  in >> _cumulativeExposureJ; 
  in >> _gType.ID; 
  if (Global::modelVersion & INCLUDES_PK_PD) {
    in >> proteomeID; 
    _proteome = ProteomeManager::getManager()->getProteome(proteomeID);
  }
  in >> boolalpha >> _SPattenuate; 
}

void DescriptiveInfection::determineDensities(int simulationTime, double cumulativeY, double ageyears, double cumulativeh, double *timeStepMaxDensity){

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
    infage=1+simulationTime-_startdate-Global::latentp;
    if ( infage >  0) {
      iduration=_duration/Global::interval;
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
          if ( Global::interval >  1) {
                normp=W_UNIFORM();
                /*
                sample the maximum density over the T-1 remaining days in the
                time interval, (where T is the duration of the time interval)
                */
                normp=pow(normp, 1.0*1/(Global::interval-1));
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
