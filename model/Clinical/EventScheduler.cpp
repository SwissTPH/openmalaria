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

#include "Clinical/EventScheduler.h"
#include "inputData.h"
#include "util/random.h"
#include "WithinHost/WithinHostModel.h"
#include "Transmission/PerHostTransmission.h"
#include "Surveys.h"
#include "util/ModelOptions.hpp"
#include "util/errors.hpp"

namespace OM { namespace Clinical {
    using namespace OM::util;

int ClinicalEventScheduler::uncomplicatedCaseDuration;
int ClinicalEventScheduler::complicatedCaseDuration;
int ClinicalEventScheduler::extraDaysAtRisk;
map<double,double> ClinicalEventScheduler::pDeathInitial;
double ClinicalEventScheduler::neg_v;


// -----  static init  -----

void ClinicalEventScheduler::init ()
{
    if (Global::interval != 1)
        throw util::xml_scenario_error ("ClinicalEventScheduler is only designed for a 1-day timestep.");
    if (! (util::ModelOptions::option (util::INCLUDES_PK_PD)))
        throw util::xml_scenario_error ("ClinicalEventScheduler requires INCLUDES_PK_PD");
    
    ESCaseManagement::init ();
}

void ClinicalEventScheduler::setParameters (const scnXml::HSEventScheduler& esData) {
    const scnXml::ClinicalOutcomes& coData = esData.getClinicalOutcomes();
    
    uncomplicatedCaseDuration = coData.getUncomplicatedCaseDuration();
    complicatedCaseDuration = coData.getComplicatedCaseDuration();
    extraDaysAtRisk = coData.getComplicatedRiskDuration() - complicatedCaseDuration;
    if( uncomplicatedCaseDuration<1 ||
	complicatedCaseDuration<1 ||
	extraDaysAtRisk+complicatedCaseDuration<1 ||
	extraDaysAtRisk>0 )
	throw util::xml_scenario_error("Clinical outcomes: constraints on case/risk duration not met (see documentation)");
    double alpha = coData.getPropDeathsFirstDay();
    if( alpha<0.0 || 1.0<alpha )
	throw util::xml_scenario_error("Clinical outcomes: Proportion of direct deaths on first day should be within range [0,1]");
    
    const map<double,double>& cfr = CaseManagementCommon::getCaseFatalityRates();
    pDeathInitial.clear();	// in case we're re-loading from intervention data
    for( map<double,double>::const_iterator it = cfr.begin(); it != cfr.end(); ++it ){
	pDeathInitial[it->first] = alpha * it->second;
    }
    neg_v = - InputData.getParameter( Params::CFR_PAR_REDUCTION_SCALAR );
}

void ClinicalEventScheduler::cleanup () {
    ESCaseManagement::cleanup ();
}


// -----  construction, destruction and checkpointing  -----

ClinicalEventScheduler::ClinicalEventScheduler (double cF, double tSF) :
        ClinicalModel (cF),
        pgState (Pathogenesis::NONE),
        timeOfRecovery (Global::TIMESTEP_NEVER)
{}
ClinicalEventScheduler::~ClinicalEventScheduler() {}


// -----  other methods  -----

void ClinicalEventScheduler::massDrugAdministration(OM::WithinHost::WithinHostModel& withinHostModel, double ageYears) {
    // Note: we use the same medication method as with drugs as treatments, hence the actual
    // medication doesn't occur until the next timestep.
    // Note: we augment any existing medications, however future medications will replace any yet-
    // to-be-medicated MDA treatments (even all MDA doses when treatment happens immediately).
    ESCaseManagement::massDrugAdministration (medicateQueue);
}

double ClinicalEventScheduler::getPDeathInitial (double ageYears) {
    assert ( ageYears >= 0.0 );
    map<double,double>::const_iterator it = pDeathInitial.upper_bound( ageYears );
    assert( it != pDeathInitial.end() );
    double a1 = it->first;
    double f1 = it->second;
    --it;
    double a0 = it->first;	// a0 <=ageYears < a1
    double f0 = it->second;
    return (ageYears - a0) / (a1 - a0) * (f1 - f0) + f0;
}

void ClinicalEventScheduler::doClinicalUpdate (WithinHost::WithinHostModel& withinHostModel, PerHostTransmission& hostTransmission, double ageYears, SurveyAgeGroup ageGroup)
{
    // Run pathogenesisModel
    // Note: we use Pathogenesis::COMPLICATED instead of Pathogenesis::SEVERE, though considering
    // the model's not designed to handle "coninfections" as presented by the Pathogenesis model
    // this is unimportant.
    Pathogenesis::State newState = pathogenesisModel->determineState (ageYears, withinHostModel);
    
    if ( Global::simulationTime == timeOfRecovery ) {
	//TODO: report sequelae?
	// if ( Sequelae )
	//     pgState = Pathogenesis::State (pgState | Pathogenesis::SEQUELAE);
	// else
	    pgState = Pathogenesis::State (pgState | Pathogenesis::RECOVERY);
	// report event, at conclusion:
	latestReport.update (Global::simulationTime, ageGroup, pgState);
	
	// Individual recovers (and is immediately susceptible to new cases)
	pgState = Pathogenesis::NONE;	// recovery (reset to healthy state)
	
	// And returns to transmission (if was removed)
	hostTransmission.removeFromTransmission( false );
    }
    
    bool cmEvent = false;	// set true when we need to do case management
    if ( pgState & Pathogenesis::COMPLICATED ) {	// previous state: complicated
	// if severe, no events happen for course of medication
    } else if ( pgState & Pathogenesis::SICK ) {	// previous state: uncomplicated (UC/UC2/NMF)
	// if uc/uc2/nmf, the only event that can happen within 3 days is progression to severe/coinfection
	if ( newState & Pathogenesis::COMPLICATED )
	    cmEvent = true;
    } else {	// previous state: healthy
	if ( newState & Pathogenesis::SICK )	// any (malarial/non-malarial) sickness
	    cmEvent = true;
    }
    
    if ( cmEvent )
    {
        if ( (newState & pgState) & Pathogenesis::MALARIA)
            newState = Pathogenesis::State (newState | Pathogenesis::SECOND_CASE);
        pgState = Pathogenesis::State (pgState | newState);
	lastEventTime = Global::simulationTime;
	
	if (pgState & Pathogenesis::MALARIA) {
	    if (util::ModelOptions::option (util::PENALISATION_EPISODES)) {
		withinHostModel.immunityPenalisation();
	    }
	}
	
	CMAuxOutput auxOut =
	    ESCaseManagement::execute (medicateQueue, pgState, withinHostModel, ageYears, ageGroup);
	
	if ( auxOut.hospitalised ) {	// in hospital
	    pgState = Pathogenesis::State (pgState | Pathogenesis::EVENT_IN_HOSPITAL);
	    
	    // Patients in hospital are removed from the transmission cycle.
	    // This should have an effect from the start of the next timestep.
	    hostTransmission.removeFromTransmission( false );
	}
	if ( auxOut.RDT_used ) {
	    Surveys.current->report_Clinical_RDTs (1);
	}
	
	if (pgState & Pathogenesis::COMPLICATED) {
	    // The probability of death on the first day is some fixed input
	    // scaled by age-specific CFR.
	    double pDeath = getPDeathInitial( ageYears );
	    if (random::uniform_01() < pDeath) {
		pgState = Pathogenesis::State (pgState | Pathogenesis::DIRECT_DEATH);
		latestReport.update (Global::simulationTime, ageGroup, pgState);
		
		// Human dies this timestep; we still allow medication of
		// today's treatments for costing.
		_doomed = DOOMED_COMPLICATED;
	    }
	    previousDensity = withinHostModel.getTotalDensity();
	    
	    // complicatedCaseDuration should to some respects be associated
	    // with medication duration, however ongoing medications after
	    // exiting hospital are OK and medications terminating before the
	    // end of hospitalization shouldn't matter too much if the person
	    // can't recieve new infections due to zero transmission in hospital.
	    timeOfRecovery = Global::simulationTime + complicatedCaseDuration;
	} else {
	    timeOfRecovery = Global::simulationTime + uncomplicatedCaseDuration;
	}
    } else {
	// No new event (haven't changed state this timestep).
	
	// Complicated case & at risk of death (note: extraDaysAtRisk <= 0)
	if( (pgState & Pathogenesis::COMPLICATED) && (Global::simulationTime < timeOfRecovery + extraDaysAtRisk) ) {
	    // In complicated episodes, S(t), the probability of survival on
	    // subsequent days t, is described by log(S(t)) = -v(Y(t)/Y(t-1)),
	    // for parasite density Y(t). v_neg below is -v.
	    double parasiteReductionEffect = withinHostModel.getTotalDensity() / previousDensity;
	    double pDeath = 1.0 - exp( neg_v * parasiteReductionEffect );
	    if (random::uniform_01() < pDeath) {
		pgState = Pathogenesis::State (pgState | Pathogenesis::DIRECT_DEATH);
		latestReport.update (Global::simulationTime, ageGroup, pgState);
		
		// Human dies this timestep; we still allow medication of
		// today's treatments for costing.
		_doomed = DOOMED_COMPLICATED;
	    }
	    previousDensity = withinHostModel.getTotalDensity();
	}
    }
    
    
    if (pgState & Pathogenesis::INDIRECT_MORTALITY && _doomed == 0)
        _doomed = -Global::interval; // start indirect mortality countdown
    
    // Process pending medications (in interal queue) and apply/update:
    for (list<MedicateData>::iterator it = medicateQueue.begin(); it != medicateQueue.end();) {
	list<MedicateData>::iterator next = it;
	++next;
	if ( it->time < 1.0 ) { // Medicate today's medications
	    withinHostModel.medicate (it->abbrev, it->qty, it->time, ageYears);
	    Surveys.current->report_Clinical_DrugUsage (it->abbrev, it->qty);
	    medicateQueue.erase (it);
	} else {   // and decrement treatment seeking delay for the rest
	    it->time -= 1.0;
	}
	it = next;
    }
}


void ClinicalEventScheduler::checkpoint (istream& stream) {
    ClinicalModel::checkpoint (stream);
    int s;
    s & stream;
    pgState = Pathogenesis::State(s);
    lastEventTime & stream;
    timeOfRecovery & stream;
    previousDensity & stream;
    medicateQueue & stream;
}
void ClinicalEventScheduler::checkpoint (ostream& stream) {
    ClinicalModel::checkpoint (stream);
    pgState & stream;
    lastEventTime & stream;
    timeOfRecovery & stream;
    previousDensity & stream;
    medicateQueue & stream;
}

} }