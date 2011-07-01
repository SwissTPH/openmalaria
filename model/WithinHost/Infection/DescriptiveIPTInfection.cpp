/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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

#include "WithinHost/Infection/DescriptiveIPTInfection.h"
#include "util/random.h"

namespace OM { namespace WithinHost {
     using namespace OM::util;
   // static IPT variables
vector<DescriptiveIPTInfection::GenotypeData> DescriptiveIPTInfection::genotypes;


// -----  static init/clear -----

// Only called if IPT is present
void DescriptiveIPTInfection::initParameters (const scnXml::IPTDescription& xmlIPTI){
  const scnXml::IPTDescription::InfGenotypeSequence& genotypesData = xmlIPTI.getInfGenotype();
  genotypes.reserve (genotypesData.size());
  
  double genotypeCumFreq = 0.0;
  for (scnXml::IPTDescription::InfGenotypeConstIterator it = genotypesData.begin(); it != genotypesData.end(); ++it) {
    genotypeCumFreq += it->getFreq();
    genotypes.push_back( GenotypeData(
        genotypeCumFreq,
        TimeStep(it->getTolPeriod()), TimeStep(it->getProph()),
        it->getACR(), it->getAtten()
    ) );
  }
  assert( genotypes.size() == genotypesData.size() );
  genotypes[genotypes.size()-1].cumFreq = 1.0;	// make sure.. for random draws
  // (otherwise we rely on specification with XML and floating-point arithmatic)
}


// -----  non-static init/destruction  -----

DescriptiveIPTInfection::DescriptiveIPTInfection(TimeStep now, TimeStep lastSPdose) :
  DescriptiveInfection(now), _SPattenuate(false)
{
    // proteome_ID is initialized to 0xFFFFFFFF
    
    double rndSample=(random::uniform_01());
    double lowerBound = 0.0;
    //This Loop assigns the infection a genotype according to its frequency
    for (size_t genotypeCounter=0; genotypeCounter < genotypes.size(); genotypeCounter++){
      if (rndSample >= lowerBound && rndSample < genotypes[genotypeCounter].cumFreq){
	proteome_ID=genotypeCounter;
	break;
      }
      lowerBound = genotypes[genotypeCounter].cumFreq;
    }
    assert (proteome_ID < genotypes.size());
    /*
    The attenuation effect of SP is only effective during a certain time-window for certain IPTi models
    If t(=now) lies within this time window, SPattenuate is true, false otherwise.
    The time window starts after the prophylactic period ended (during the prophylactic
    period infections are cleared) and ends genotypeTolPeriod(iTemp%iData%gType%ID) time steps later.
    */
    if (now-lastSPdose > genotypes[proteome_ID].proph &&
	now-lastSPdose <= genotypes[proteome_ID].proph + genotypes[proteome_ID].tolPeriod){
      _SPattenuate=true;
    }
}

bool DescriptiveIPTInfection::eventSPClears (TimeStep _lastSPDose) {
    if(TimeStep::simulation - _startdate < latentp)
	return false;	// don't consider pre-patent infections
    
    return
	(TimeStep::simulation - _lastSPDose <= DescriptiveIPTInfection::genotypes[proteome_ID].proph)
	&& (random::uniform_01() <= DescriptiveIPTInfection::genotypes[proteome_ID].ACR);
}

double DescriptiveIPTInfection::asexualAttenuation () {
  double attFact = 1.0 / genotypes[proteome_ID].atten;
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