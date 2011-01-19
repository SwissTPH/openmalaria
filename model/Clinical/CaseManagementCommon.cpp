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
#include "Clinical/ImmediateOutcomes.h"
#include "Clinical/ESCaseManagement.h"
#include "inputData.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
#include <limits>
#include <boost/format.hpp>

namespace OM { namespace Clinical {
    using util::AgeGroupInterpolation;
    
    TimeStep CaseManagementCommon::healthSystemSource( TimeStep::never );
    util::AgeGroupInterpolation* CaseManagementCommon::caseFatalityRate = AgeGroupInterpolation::dummyObject();
    double CaseManagementCommon::_oddsRatioThreshold;
    util::AgeGroupInterpolation* CaseManagementCommon::pSeqInpatient = AgeGroupInterpolation::dummyObject();
    
    // -----  functions  -----
    
    void CaseManagementCommon::initCommon (){
	_oddsRatioThreshold = exp (InputData.getParameter (Params::LOG_ODDS_RATIO_CF_COMMUNITY));
	
	changeHealthSystem(TimeStep::never);
    }
    void CaseManagementCommon::cleanupCommon (){
        AgeGroupInterpolation::freeObject( caseFatalityRate );
        AgeGroupInterpolation::freeObject( pSeqInpatient );
    }
    
    void CaseManagementCommon::changeHealthSystem (TimeStep source) {
	healthSystemSource = source;
	const scnXml::HealthSystem& healthSystem = getHealthSystem ();
	
	readCommon (healthSystem);
	
	if (util::ModelOptions::option (util::CLINICAL_EVENT_SCHEDULER))
	    ESCaseManagement::setHealthSystem(healthSystem);
	else
	    ClinicalImmediateOutcomes::setHealthSystem(healthSystem);
    }
    
    const scnXml::HealthSystem& CaseManagementCommon::getHealthSystem () {
	if (healthSystemSource == TimeStep::never) {
	    return InputData().getHealthSystem();
	} else {
	    const scnXml::Intervention* interv = InputData.getInterventionByTime (healthSystemSource);
	    if (interv == NULL || !interv->getChangeHS().present())
		throw runtime_error ("healthSystemSource invalid");
	    return interv->getChangeHS().get();
	}
	assert(false);	// unreachable
    }
    
    void CaseManagementCommon::readCommon (const scnXml::HealthSystem& healthSystem)
    {
	// -----  case fatality rates  -----
        AgeGroupInterpolation::freeObject( caseFatalityRate );
	caseFatalityRate = AgeGroupInterpolation::makeObject( healthSystem.getCFR(), "CFR" );
	
	// -----  sequelae  -----
        AgeGroupInterpolation::freeObject( pSeqInpatient );
        pSeqInpatient = AgeGroupInterpolation::makeObject( healthSystem.getPSequelaeInpatient(), "pSequelaeInpatient" );
    }
    
    double CaseManagementCommon::getCommunityCaseFatalityRate (double caseFatalityRatio)
    {
	double x = caseFatalityRatio * _oddsRatioThreshold;
	return x / (1 - caseFatalityRatio + x);
    }

    void CaseManagementCommon::staticCheckpoint (ostream& stream) {
	healthSystemSource & stream;
    }
    void CaseManagementCommon::staticCheckpoint (istream& stream) {
	healthSystemSource & stream;
	if (healthSystemSource != TimeStep::never)
	    changeHealthSystem( healthSystemSource );
    }
} }
