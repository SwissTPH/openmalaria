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

#include "Clinical/CaseManagementCommon.h"
#include "Clinical/OldCaseManagement.h"
#include "Clinical/ESCaseManagement.h"
#include "inputData.h"
#include "util/ModelOptions.hpp"
#include "util/errors.hpp"
#include <limits>
#include <boost/format.hpp>

namespace OM { namespace Clinical {
    int CaseManagementCommon::healthSystemSource;
    map<double,double> CaseManagementCommon::caseFatalityRates;
    
    // -----  functions  -----
    
    void CaseManagementCommon::changeHealthSystem (int source) {
	healthSystemSource = source;
	const scnXml::HealthSystem& healthSystem = getHealthSystem ();
	
	if (util::ModelOptions::option (util::CLINICAL_EVENT_SCHEDULER))
	    ESCaseManagement::setHealthSystem(healthSystem);
	else
	    OldCaseManagement::setHealthSystem(healthSystem);
	
	readCaseFatalityRatio (healthSystem);
    }
    
    const scnXml::HealthSystem& CaseManagementCommon::getHealthSystem () {
	if (healthSystemSource == -1) {
	    return InputData().getHealthSystem();
	} else {
	    const scnXml::Intervention* interv = InputData.getInterventionByTime (healthSystemSource);
	    if (interv == NULL || !interv->getChangeHS().present())
		throw runtime_error ("healthSystemSource invalid");
	    return interv->getChangeHS().get();
	}
	assert(false);	// unreachable
    }
    
    void CaseManagementCommon::readCaseFatalityRatio (const scnXml::HealthSystem& healthSystem)
    {
	caseFatalityRates.clear();	// Necessary when re-read from an intervention
	
	const scnXml::CFRAgeGroups::GroupSequence& xmlGroupSeq = healthSystem.getCFR().getGroup();
	
	BOOST_FOREACH( const scnXml::Group& xmlGroup, xmlGroupSeq ){
	    double lbound = xmlGroup.getLowerbound();
	    if( !caseFatalityRates.insert( make_pair( lbound, xmlGroup.getCfr() ) ).second )
		throw util::xml_scenario_error( (boost::format("CFR: lower bound %1% listed twice") %lbound).str() );
	}
	// CFR is constant for everyone above the highest non-inf upperbound
	caseFatalityRates[ numeric_limits<double>::infinity() ] =
	    caseFatalityRates.rbegin()->second;
	
	// first lower-bound must be 0
	if( caseFatalityRates.begin()->first != 0.0 )
	    throw util::xml_scenario_error( "CFR: first lower-bound must be 0" );
    }
    
    double CaseManagementCommon::caseFatality (double ageYears)
    {
	assert ( ageYears >= 0.0 );
	map<double,double>::const_iterator it = caseFatalityRates.upper_bound( ageYears );
	assert( it != caseFatalityRates.end() );
	double a1 = it->first;
	double f1 = it->second;
	--it;
	double a0 = it->first;	// a0 <=ageYears < a1
	double f0 = it->second;
	return (ageYears - a0) / (a1 - a0) * (f1 - f0) + f0;
    }
    
    void CaseManagementCommon::staticCheckpoint (ostream& stream) {
	healthSystemSource & stream;
    }
    void CaseManagementCommon::staticCheckpoint (istream& stream) {
	healthSystemSource & stream;
	if (healthSystemSource != -1)
	    changeHealthSystem( healthSystemSource );
    }
} }