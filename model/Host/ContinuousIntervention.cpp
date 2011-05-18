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
#include "Host/Human.h"
#include "inputData.h"
#include "util/random.h"
#include "util/errors.h"
#include <algorithm>
#include <cmath>

namespace OM { namespace Host {

vector<ContinuousIntervention::AgeIntervention> ContinuousIntervention::ctsIntervs;

void ContinuousIntervention::init (
	void (Human::*deployVaccine) (),
	void (Human::*deployItn) (),
	void (Human::*deployIpti) (),
	void (Human::*deployCohort) ()
) {
    /* FIXME
    // NOTE: "vaccine" interventions we don't deal with; they have some rather
    // special deployment rules (see Human::_lastVaccineDose).
    
    const scnXml::ContinuousInterv xmlCts = InputData().getInterventions().getContinuous().get();
    const scnXml::ContinuousInterv::VaccineSequence& seqVaccine = xmlCts.getVaccine();
    const scnXml::ContinuousInterv::ITNSequence& seqItn = xmlCts.getITN();
    const scnXml::ContinuousInterv::IptiSequence& seqIpti = xmlCts.getIpti();
    const scnXml::ContinuousInterv::CohortSequence& seqCohort = xmlCts.getCohort();
    size_t n = seqVaccine.size() + seqItn.size() + seqIpti.size() + seqCohort.size();
    ctsIntervs.reserve( n );
    
    typedef xsd::cxx::tree::sequence< ::scnXml::AgeSpecific >::const_iterator const_it_t;
    
    // Now read each type of intervention in turn, then sort.
    // Changing XML structure could reduce code duplication here.
    
    for (const_it_t it = seqVaccine.begin(); it != seqVaccine.end(); ++it){
	ctsIntervs.push_back(AgeIntervention( *it, deployVaccine ));
    }
    for (const_it_t it = seqItn.begin(); it != seqItn.end(); ++it){
        ctsIntervs.push_back(AgeIntervention( *it, deployItn ));
    }
    for (const_it_t it = seqIpti.begin(); it != seqIpti.end(); ++it){
        ctsIntervs.push_back(AgeIntervention( *it, deployIpti ));
    }
    for (const_it_t it = seqCohort.begin(); it != seqCohort.end(); ++it){
        ctsIntervs.push_back(AgeIntervention( *it, deployCohort ));
    }
    
    sort( ctsIntervs.begin(), ctsIntervs.end() );
    assert( ctsIntervs.size() == n );
    */
}

ContinuousIntervention::AgeIntervention::AgeIntervention(
    const ::scnXml::AgeSpecific& elt, void(Human::*func) ()
) :
    begin( elt.getBegin() ),
    end( elt.getEnd() ),
    ageTimesteps( TimeStep::fromYears( elt.getTargetAgeYrs() ) ),
    cohortOnly( elt.getCohort() ),
    coverage( elt.getCoverage() ),
    deploy( func )
{
    if( ageTimesteps <= TimeStep(0) ){
	ostringstream msg;
	msg << "continuous intervention with target age "<<elt.getTargetAgeYrs();
	msg << " years corresponds to timestep "<<ageTimesteps;
	msg << "; must be at least timestep 1.";
	throw util::xml_scenario_error( msg.str() );
    }
}

void ContinuousIntervention::deploy (Human* human, TimeStep ageTimesteps){
    while( nextCtsDist < ctsIntervs.size() ){
	if( ctsIntervs[nextCtsDist].ageTimesteps > ageTimesteps )
	    break;	// remaining intervs happen in future
	// If interv for now, do it. (If we missed the time, ignore it.)
	if( ctsIntervs[nextCtsDist].ageTimesteps == ageTimesteps ){
	    if( ctsIntervs[nextCtsDist].begin <= TimeStep::interventionPeriod && TimeStep::interventionPeriod <= ctsIntervs[nextCtsDist].end ){
		if( !ctsIntervs[nextCtsDist].cohortOnly || human->getInCohort() ){
		    if (util::random::uniform_01() < ctsIntervs[nextCtsDist].coverage){
			(human->*(ctsIntervs[nextCtsDist].deploy)) ();
		    }
		}
	    }
	}
	nextCtsDist++;
    }
}

} }
