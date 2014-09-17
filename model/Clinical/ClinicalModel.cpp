/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2014 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2014 Liverpool School Of Tropical Medicine
 * 
 * OpenMalaria is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "Clinical/ClinicalModel.h"

#include "Clinical/CaseManagementCommon.h"
#include "Clinical/EventScheduler.h"
#include "Clinical/ImmediateOutcomes.h"
#include "Clinical/DecisionTree5Day.h"
#include "Host/NeonatalMortality.h"

#include "Monitoring/Survey.h"
#include "util/ModelOptions.h"
#include <schema/scenario.h>

namespace OM { namespace Clinical {
    using namespace Monitoring;

bool opt_event_scheduler = false;
bool opt_imm_outcomes = false;

// -----  static methods  -----

void ClinicalModel::init( const Parameters& parameters, const scnXml::Scenario& scenario ) {
    const scnXml::Clinical& clinical = scenario.getModel().getClinical();
    //FIXME(schema): input should be in days
    SimTime hsMemory = sim::fromTS(TimeStep(clinical.getHealthSystemMemory()));
    initCMCommon( parameters, hsMemory );
    
    if (util::ModelOptions::option (util::CLINICAL_EVENT_SCHEDULER)){
        opt_event_scheduler = true;
        ClinicalEventScheduler::init( parameters, clinical );
    }else{
        if( scenario.getHealthSystem().getImmediateOutcomes().present() )
            opt_imm_outcomes = true;
        // else: decision tree 5 day
    }
}

void ClinicalModel::changeHS( const scnXml::HealthSystem& healthSystem ){
    caseFatalityRate.set( healthSystem.getCFR(), "CFR" );
    pSequelaeInpatient.set( healthSystem.getPSequelaeInpatient(), "pSequelaeInpatient" );
    if( opt_event_scheduler ){
        if( !healthSystem.getEventScheduler().present() ){
            throw util::xml_scenario_error ("Expected EventScheduler "
                "section in healthSystem data (initial or intervention)");
        }
        ESCaseManagement::setHealthSystem(healthSystem.getEventScheduler().get());
    }else if( opt_imm_outcomes ){
        if ( !healthSystem.getImmediateOutcomes().present() ){
            throw util::xml_scenario_error ("Expected ImmediateOutcomes "
                "section in healthSystem data (initial or intervention)");
        }
        ImmediateOutcomes::setHealthSystem(healthSystem.getImmediateOutcomes().get());
    }else{
        if( !healthSystem.getDecisionTree5Day().present() ){
            throw util::xml_scenario_error ("Expected DecisionTree5Day "
                "section in healthSystem data (initial or intervention)");
        }
        DecisionTree5Day::setHealthSystem(healthSystem.getDecisionTree5Day().get());
    }
}

ClinicalModel* ClinicalModel::createClinicalModel (double tSF) {
    if (opt_event_scheduler){
        return new ClinicalEventScheduler (tSF);
    }else if( opt_imm_outcomes ){
        return new ImmediateOutcomes (tSF);
    }else{
        return new DecisionTree5Day( tSF );
    }
}


// -----  non-static construction, destruction and checkpointing  -----

ClinicalModel::ClinicalModel () :
    doomed(NOT_DOOMED)
{}
ClinicalModel::~ClinicalModel () {
  // latestReport is reported, if any, by destructor
}


// -----  other non-static methods  -----

bool ClinicalModel::isDead (TimeStep ageTimeSteps) {
    if (ageTimeSteps > TimeStep::maxAgeIntervals)	// too old (reached age limit)
        doomed = DOOMED_TOO_OLD;
    if (doomed > NOT_DOOMED)	// killed by some means
        return true;	// remove from population
    return false;
}

void ClinicalModel::update (Human& human, double ageYears, bool newBorn) {
    if (doomed < NOT_DOOMED)	// Countdown to indirect mortality
        doomed -= TimeStep::interval;
    
    //indirect death: if this human's about to die, don't worry about further episodes:
    if (doomed <= DOOMED_EXPIRED) {	//clinical bout 6 intervals before
        Survey::current().addInt( Report::MI_INDIRECT_DEATHS, human, 1 );
        doomed = DOOMED_INDIRECT;
        return;
    }
    if(newBorn /* i.e. first update since birth */) {
        // Chance of neonatal mortality:
        if (Host::NeonatalMortality::eventNeonatalMortality()) {
            Survey::current().addInt( Report::MI_INDIRECT_DEATHS, human, 1 );
            doomed = DOOMED_NEONATAL;
            return;
        }
    }
    
    doClinicalUpdate (human, ageYears);
}

void ClinicalModel::updateInfantDeaths( TimeStep ageTimeSteps ){
    // update array for the infant death rates
    if (ageTimeSteps <= TimeStep::intervalsPerYear){
        infantIntervalsAtRisk[ageTimeSteps.asInt()-1] += 1;     // baseline
        // Testing doomed == DOOMED_NEXT_TS gives very slightly different results than
        // testing doomed == DOOMED_INDIRECT (due to above if(..))
        if( doomed == DOOMED_COMPLICATED || doomed == DOOMED_NEXT_TS || doomed == DOOMED_NEONATAL ){
            infantDeaths[ageTimeSteps.asInt()-1] += 1;  // deaths
        }
    }
}


void ClinicalModel::checkpoint (istream& stream) {
    latestReport & stream;
    doomed & stream;
}
void ClinicalModel::checkpoint (ostream& stream) {
    latestReport & stream;
    doomed & stream;
}

} }