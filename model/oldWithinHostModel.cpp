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

#include "GSLWrapper.h"
#include "human.h"
#include "Infection.h"
#include "oldWithinHostModel.h"
#include "simulation.h"
#include "intervention.h"

using namespace std;

OldWithinHostModel::OldWithinHostModel(Human *human) : WithinHostModel(human) {
  _SPattenuationt = 0;
}

void OldWithinHostModel::calculateDensity(Infection *inf) {
  //Age of infection. (Blood stage infection starts latentp intervals later than inoculation ?)
  int infage=1+Simulation::simulationTime-inf->getStartDate()-Global::latentp;
  
  if ( infage >  0) {
    double y;
    if ( infage <=  maxDur) {
      int iduration=inf->getDuration()/Global::interval;
      if ( iduration >  maxDur)
        iduration=maxDur;
      
      y=(float)exp(inf->getMeanLogParasiteCount(infage - 1 + (iduration - 1)*maxDur));
    } else {
      y=(float)exp(inf->getMeanLogParasiteCount(maxDur - 1 + (maxDur - 1)*maxDur));
    }
    if ( y <  1.0)
      y=1.0;
    
    //effect of cumulative Parasite density (named Dy in AJTM)
    double dY;
    //effect of number of infections experienced since birth (named Dh in AJTM)
    double dH;
    //effect of age-dependent maternal immunity (named Dm in AJTM)
    double dA;
    
    if ( cumulativeh <=  1.0) {
      dY=1;
      dH=1;
    } else {
      dH=1 / (1+(cumulativeh-1.0)/inf->getCumulativeHstar());
      //TODO: compare this with the asex paper
      dY=1 / (1+(cumulativeY-inf->getCumulativeExposureJ())/inf->getCumulativeYstar());
    }
    dA=1 - inf->getAlpha_m() * exp(-inf->getDecayM() * _human->getAgeInYears());
    
    double survival=min (dY*dH*dA, 1.0);
    double logy=log(y)*survival;
    /*
      The expected parasite density in the non naive host. 
      As regards the second term in AJTM p.9 eq. 9, in published and current implementations Dx is zero.
    */
    y=exp(logy);
    //Perturb y using a lognormal 
    double varlog=inf->getSigma0sq()/(1+(cumulativeh/inf->getXNuStar()));
    double stdlog=sqrt(varlog);
    /*
      This code samples from a log normal distribution with mean equal to the predicted density
      n.b. AJTM p.9 eq 9 implies that we sample the log of the density from a normal with mean equal to
      the log of the predicted density.  If we really did the latter then this bias correction is not needed.
    */
    double meanlog=log(y)-stdlog*stdlog/2.0;
    timeStepMaxDensity = 0.0;
    if ( stdlog >  0.0000001) {
      if ( Global::interval >  1) {
	double normp=W_UNIFORM();
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
	timeStepMaxDensity = sampleFromLogNormal(normp, meanlog, stdlog);
      }
      //calculate the expected density on the day of sampling
      y=(float)sampleFromLogNormal(W_UNIFORM(), meanlog, stdlog);
      timeStepMaxDensity = std::max( y, timeStepMaxDensity);
    }
    if (( y >  maxDens) || ( timeStepMaxDensity >  (double)maxDens)) {
      cout << "MD lim" << endl;
      y=maxDens;
      timeStepMaxDensity = (double) y;
    }
    inf->setDensity(y);
  }
  else {
    inf->setDensity(0.0);
  }
}

void OldWithinHostModel::calculateDensities() {
  //double ageyears = _human->getAgeInYears();
  _human->setPTransmit(0);
  _human->setPatentInfections(0);
  _human->setTotalDensity(0.0);
  _human->setTimeStepMaxDensity(0.0);
  if ( _human->getCumulativeInfections() >  0) {
    cumulativeh=_human->getCumulativeh();
    cumulativeY=_human->getCumulativeY();
    // IPTi SP dosec lears infections at the time that blood-stage parasites appear     
    if (IPTIntervention::IPT) {
      SPAction();
    }
    std::list<Infection*>::iterator i;
    for(i=_human->getInfections()->begin(); i!=_human->getInfections()->end(); i++){
      //std::cout<<"uis: "<<infData->duration<<std::endl;
      timeStepMaxDensity=_human->getTimeStepMaxDensity();
      if (Global::modelVersion & WITHIN_HOST_PARASITE) {
        (*i)->setDensity((*i)->determineWithinHostDensity());
        timeStepMaxDensity=std::max((double)(*i)->getDensity(), timeStepMaxDensity);
        _human->setTimeStepMaxDensity(timeStepMaxDensity);
      }
      else {
        if (Global::modelVersion & MAX_DENS_RESET) {
          timeStepMaxDensity=0.0;
        }
	calculateDensity(*i);
        //(*i)->determineDensities(Simulation::simulationTime, _human->getCumulativeY(), ageyears, cumulativeh , &(timeStepMaxDensity));
        (*i)->setDensity((double)(*i)->getDensity()*exp(-_human->getInnateImmunity()));

        /*
          Possibly a better model version ensuring that the effect of variation in innate immunity
          is reflected in case incidence would have the following here:
        */
        if (Global::modelVersion & INNATE_MAX_DENS) {
          timeStepMaxDensity=(double)timeStepMaxDensity*exp(-_human->getInnateImmunity());
        }
        //Include here the effect of blood stage vaccination
        if (Vaccine::BSV.active) {
          (*i)->setDensity((*i)->getDensity()*(1-_human->getBSVEfficacy()));
          timeStepMaxDensity=(double)timeStepMaxDensity*(1-_human->getBSVEfficacy());
        }
        // Include here the effect of attenuated infections by SP concentrations
        if (Global::modelVersion & ATTENUATION_ASEXUAL_DENSITY) {
          if ( IPTIntervention::IPT &&  (*i)->getSPattenuate() ==  1) {
            (*i)->setDensity((*i)->getDensity()/IPTIntervention::genotypeAtten[(*i)->getGenoTypeID() - 1]);
            timeStepMaxDensity=(double)timeStepMaxDensity/IPTIntervention::genotypeAtten[(*i)->getGenoTypeID() - 1];
            _SPattenuationt=(int)std::max(_SPattenuationt*1.0, ((*i)->getStartDate()+((*i)->getDuration()/Global::interval) * IPTIntervention::genotypeAtten[(*i)->getGenoTypeID() - 1]));
          }
        }
        if (Global::modelVersion & MAX_DENS_CORRECTION) {
          _human->setTimeStepMaxDensity(std::max(timeStepMaxDensity, _human->getTimeStepMaxDensity()));
        }
        else {
          _human->setTimeStepMaxDensity(timeStepMaxDensity);
        }
      }
      _human->setTotalDensity(_human->getTotalDensity()+(*i)->getDensity());
      //Compute the proportion of parasites remaining after innate blood stage effect
      if ( (*i)->getDensity() > Human::detectionlimit) {
        _human->setPatentInfections(_human->getPatentInfections()+1);
      }
      if ( (*i)->getStartDate() == (Simulation::simulationTime-1)) {
        _human->setCumulativeh(_human->getCumulativeh()+1);
      }
      (*i)->setDensity(std::min(maxDens, (*i)->getDensity()));
      (*i)->setCumulativeExposureJ((*i)->getCumulativeExposureJ()+Global::interval*(*i)->getDensity());
      _human->setCumulativeY(_human->getCumulativeY()+Global::interval*(*i)->getDensity());
    }
    if (Global::modelVersion & ATTENUATION_ASEXUAL_DENSITY) {
      if ( IPTIntervention::IPT &&  _SPattenuationt > Simulation::simulationTime &&  _human->getTotalDensity() <  10) {
        _human->setTotalDensity(10);
        _human->setCumulativeY(_human->getCumulativeY()+10);
      }
    }
  }
  _human->setPTransmit(_human->infectiousness());
}

void OldWithinHostModel::SPAction(){
 
  /*TODO if we want to look at presumptive SP treatment with the PkPD model we
    need to add some code here that will be conditionally implemented depending on the
    model version.*/

  double rnum;
  std::list<Infection*>::iterator i=_human->getInfections()->begin();
  while(i != _human->getInfections()->end()){
    if ( 1+Simulation::simulationTime-(*i)->getStartDate()-Global::latentp > 0){
      rnum=W_UNIFORM();
      if ((rnum<=IPTIntervention::genotypeACR[(*i)->getGenoTypeID()-1]) &&
           (Simulation::simulationTime - _human->getLastSPDose() <= IPTIntervention::genotypeProph[(*i)->getGenoTypeID()-1])) {
        delete *i;
        i=_human->getInfections()->erase(i);
        _human->setMOI(_human->getMOI()-1);
      }
      else{
        i++;
      }
    }
    else{
      i++;
    }
  }
}

void OldWithinHostModel::read(istream& in) {
  in >> _SPattenuationt;
  in >> cumulativeY;
  in >> cumulativeh;
  in >> timeStepMaxDensity;
}

void OldWithinHostModel::write(ostream& out) const {
  out << _SPattenuationt << endl;
  out << cumulativeY << endl;
  out << cumulativeh << endl;
  out << timeStepMaxDensity << endl;
}
