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

#include "Host/ContinuousIntervention.h"
#include "inputData.h"
#include "util/random.h"
#include <algorithm>
#include <cmath>

namespace OM { namespace Host {

vector<ContinuousIntervention::AgeIntervention> ContinuousIntervention::ctsIntervs;

void ContinuousIntervention::init (
	void (Human::*deployItn) (),
	void (Human::*deployIpti) ()
) {
    if (!InputData().getInterventions().getContinuous().present()) {
	return;
    }
    
    // NOTE: "vaccine" interventions we don't deal with; they have some rather
    // special deployment rules (see Human::_lastVaccineDose).
    
    const scnXml::Continuous xmlCts = InputData().getInterventions().getContinuous().get();
    const scnXml::Continuous::ITNSequence& seqItn = xmlCts.getITN();
    const scnXml::Continuous::IptiSequence& seqIpti = xmlCts.getIpti();
    size_t n = seqItn.size() + seqIpti.size();
    ctsIntervs.resize( n );
    
    // Now read each type of intervention in turn, then sort.
    // Changing XML structure could reduce code duplication here.
    n = 0;	// use as index in ctsIntervs
    
    for (scnXml::Continuous::ITNConstIterator it = seqItn.begin(); it != seqItn.end(); ++it){
	ctsIntervs[n].ageTimesteps = static_cast<uint32_t>(
	    floor( it->getTargetAgeYrs() * Global::DAYS_IN_YEAR / (1.0*Global::interval) )
	);
	ctsIntervs[n].coverage = it->getCoverage();
	ctsIntervs[n].deploy = deployItn;
	n++;
    }
    for (scnXml::Continuous::IptiConstIterator it = seqIpti.begin(); it != seqIpti.end(); ++it){
	ctsIntervs[n].ageTimesteps = static_cast<uint32_t>(
	    floor( it->getTargetAgeYrs() * Global::DAYS_IN_YEAR / (1.0*Global::interval) )
	);
	ctsIntervs[n].coverage = it->getCoverage();
	ctsIntervs[n].deploy = deployIpti;
	n++;
    }
    
    sort( ctsIntervs.begin(), ctsIntervs.end() );
}

void ContinuousIntervention::deploy (Human* human, int ageTimesteps){
    while( nextCtsDist < ctsIntervs.size() ){
	if( ctsIntervs[nextCtsDist].ageTimesteps > ageTimesteps )
	    break;	// remaining intervs happen in future
	// If interv for now, do it. (If we missed the time, ignore it.)
	if( ctsIntervs[nextCtsDist].ageTimesteps == ageTimesteps ){
	    if (util::random::uniform_01() < ctsIntervs[nextCtsDist].coverage){
		(human->*(ctsIntervs[nextCtsDist].deploy)) ();
	    }
	}
	nextCtsDist++;
    }
}

} }
