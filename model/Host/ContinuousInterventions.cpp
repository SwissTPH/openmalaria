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

#include "Host/ContinuousInterventions.h"
#include "inputData.h"

vector<InterventionsPerAge> ContinuousInterventions::intervs;


//InterventionsPerAge::

void ContinuousInterventions::deploy (int ageTS) {
    if (next < intervs.size() && intervs[next].getAgeTimeSteps() == ageTS) {
	
	//...
	
	++next;
    }
}


void ContinuousInterventions::initParameters () {
    const scnXml::Interventions::ContinuousOptional& ctsXml = getInterventions().getContinuous();
    if (!ctsXml.present()) return;
    
    map<int,InterventionsPerAge> parsedIntervs;
    
    const scnXml::Continuous::VaccineSequence& vseq = ctsXml.get().getVaccine();
    for (scnXml::Continuous::VaccineConstIterator itv = vseq.begin(); itv != vseq.end(); ++itv) {
	int ageTS = (int) floor (itv->getTargetAgeYrs() * daysInYear / (double) Global::interval);
	if (parsedIntervs.count(ageTS) == 0)
	    parsedIntervs[ageTS] = InterventionsPerAge(ageTS);
	parsedIntervs[ageTS].addVaccine (itv->getCoverage());
    }
    
    const scnXml::Continuous::ITNSequence& iseq = ctsXml.get().getITN();
    for (scnXml::Continuous::ITNConstIterator itv = iseq.begin(); itv != iseq.end(); ++itv) {
	int ageTS = (int) floor (itv->getTargetAgeYrs() * daysInYear / (double) Global::interval);
	if (parsedIntervs.count(ageTS) == 0)
	    parsedIntervs[ageTS] = InterventionsPerAge(ageTS);
	parsedIntervs[ageTS].addITN (itv->getCoverage());
    }
    
    const scnXml::Continuous::IptiSequence& ipseq = ctsXml.get().getVaccine();
    for (scnXml::Continuous::IptiConstIterator itv = ipseq.begin(); itv != ipseq.end(); ++itv) {
	int ageTS = (int) floor (itv->getTargetAgeYrs() * daysInYear / (double) Global::interval);
	if (parsedIntervs.count(ageTS) == 0)
	    parsedIntervs[ageTS] = InterventionsPerAge(ageTS);
	parsedIntervs[ageTS].addIPTI (itv->getCoverage());
    }
    
    intervs.reserve(parsedIntervs.size());
    for (map<int,InterventionsPerAge>::const_iterator it = parsedIntervs.begin(); it != parsedIntervs.end(); ++it)
	intervs.push_back (it->second);
}
