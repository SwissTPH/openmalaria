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

#include "WithinHost/DescriptiveIPTInfection.h"
#include "inputData.h"
#include "util/gsl.h"

namespace OM { namespace WithinHost {
    // static IPT variables
int DescriptiveIPTInfection::numberOfGenoTypes;
double *DescriptiveIPTInfection::genotypeFreq;
int *DescriptiveIPTInfection::genotypeProph;
int *DescriptiveIPTInfection::genotypeTolPeriod;
double *DescriptiveIPTInfection::genotypeACR;
double *DescriptiveIPTInfection::genotypeAtten;


// -----  static init/clear -----

// Only called if IPT is present
void DescriptiveIPTInfection::initParameters (const scnXml::Interventions& xmlInterventions){
  const scnXml::IptDescription& xmlIPTI = xmlInterventions.getIptiDescription().get();
  
  const scnXml::IptDescription::InfGenotypeSequence& genotypes = xmlIPTI.getInfGenotype();
  numberOfGenoTypes = genotypes.size();
  genotypeFreq	= (double*)malloc(((numberOfGenoTypes))*sizeof(double));
  genotypeACR	= (double*)malloc(((numberOfGenoTypes))*sizeof(double));
  genotypeProph	= (int*)malloc(((numberOfGenoTypes))*sizeof(int));
  genotypeTolPeriod = (int*)malloc(((numberOfGenoTypes))*sizeof(int));
  genotypeAtten	= (double*)malloc(((numberOfGenoTypes))*sizeof(double));
  
  size_t i = 0;
  for (scnXml::IptDescription::InfGenotypeConstIterator it = genotypes.begin();
       it != genotypes.end(); ++it, ++i) {
    genotypeFreq[i]	= it->getFreq();
    genotypeACR[i]	= it->getACR();
    genotypeProph[i]	= it->getProph();
    genotypeTolPeriod[i]= it->getTolPeriod();
    genotypeAtten[i]	= it->getAtten();
  }
}

void DescriptiveIPTInfection::clearParameters () {
  free(genotypeFreq);
  free(genotypeACR);
  free(genotypeProph);
  free(genotypeTolPeriod);
  free(genotypeAtten);
}


// -----  non-static init/destruction  -----

DescriptiveIPTInfection::DescriptiveIPTInfection(int lastSPdose) :
  DescriptiveInfection(), _SPattenuate(false)
{
    // proteome_ID is initialized to 0xFFFFFFFF
  
  double uniformRandomVariable=(gsl::rngUniform());
    double lowerIntervalBound=0.0;
    double upperIntervalBound=genotypeFreq[0];
    //This Loop assigns the infection a genotype according to its frequency
    for (int genotypeCounter=1; genotypeCounter<=numberOfGenoTypes; genotypeCounter++){
	//FIXME: use >= and < (make sure all potential values of uniformRandomVariable are covered)
      if (uniformRandomVariable > lowerIntervalBound &&
	uniformRandomVariable < upperIntervalBound){
	proteome_ID=genotypeCounter - 1;
      }
      lowerIntervalBound=upperIntervalBound;
      /*
      The loop should never come into this else-part (exits before), so the if statement is not necessary.
      For safety reason we do it nevertheless
      */
      if ( genotypeCounter !=  numberOfGenoTypes) {
	  //FIXME: array access is out of bounds when genotypeCounter == numberOfGenoTypes
	upperIntervalBound = upperIntervalBound + genotypeFreq[genotypeCounter];
      }
      else {
	upperIntervalBound=1.0;
	//FIXME: this check wasn't here before. I've just found it fails. DH
	//assert (false);	// 'should not be here'
      }
    }
    assert (proteome_ID < 0xFFFFFFFF);
    /*
    The attenuation effect of SP is only effective during a certain time-window for certain IPTi models
    If t(=now) lies within this time window, SPattenuate is true, false otherwise.
    The time window starts after the prophylactic period ended (during the prophylactic
    period infections are cleared) and ends genotypeTolPeriod(iTemp%iData%gType%ID) time steps later.
    */
    if (Global::simulationTime-lastSPdose > genotypeProph[proteome_ID] &&
	Global::simulationTime-lastSPdose <= genotypeProph[proteome_ID] + genotypeTolPeriod[proteome_ID]){
      _SPattenuate=true;
    }
}


double DescriptiveIPTInfection::asexualAttenuation () {
  double attFact = 1.0 / genotypeAtten[proteome_ID];
  _density *= attFact;
  return attFact;
}


DescriptiveIPTInfection::DescriptiveIPTInfection (istream& stream) :
    DescriptiveInfection(stream)
{
  _SPattenuate & stream; 
}
void DescriptiveIPTInfection::checkpoint (ostream& stream) {
    DescriptiveInfection::checkpoint (stream);
    _SPattenuate & stream; 
}

} }