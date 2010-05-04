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

#include <limits>

namespace OM { namespace Clinical {
    using namespace OM::util;

int ClinicalEventScheduler::maxUCSeekingMemory;
int ClinicalEventScheduler::uncomplicatedCaseDuration;
int ClinicalEventScheduler::complicatedCaseDuration;
int ClinicalEventScheduler::extraDaysAtRisk;
double ClinicalEventScheduler::pImmediateUC;
map<double,double> ClinicalEventScheduler::pDeathInitial;
double ClinicalEventScheduler::neg_v;
double ClinicalEventScheduler::communityOddsMultiplier;


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
    
    maxUCSeekingMemory = coData.getMaxUCSeekingMemory();
    uncomplicatedCaseDuration = coData.getUncomplicatedCaseDuration();
    complicatedCaseDuration = coData.getComplicatedCaseDuration();
    extraDaysAtRisk = coData.getComplicatedRiskDuration() - complicatedCaseDuration;
    if( uncomplicatedCaseDuration<1
	|| complicatedCaseDuration<1
	|| maxUCSeekingMemory<0
	|| extraDaysAtRisk+complicatedCaseDuration<1	// at risk at least 1 day
	|| extraDaysAtRisk>0	// at risk longer than case duration
    )
	throw util::xml_scenario_error("Clinical outcomes: constraints on case/risk/memory duration not met (see documentation)");
    pImmediateUC = coData.getPImmediateUC();
    double alpha = coData.getPropDeathsFirstDay();
    communityOddsMultiplier = coData.getCommunityOddsMultiplier();
    if( !(0.0<=alpha && alpha<=1.0)
	|| !(0.0<=pImmediateUC && pImmediateUC<=1.0)
    )
	throw util::xml_scenario_error("Clinical outcomes: pImmediateUC and propDeathsFirstDay should be within range [0,1]");
    if( !(communityOddsMultiplier >= 0) )
	throw util::xml_scenario_error("Clinical outcomes: communityOddsMultiplier must be >= 0");
    
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
        caseStartTime (Global::TIMESTEP_NEVER),
        timeOfRecovery (Global::TIMESTEP_NEVER),
        timeLastTreatment (Global::TIMESTEP_NEVER),
        previousDensity (numeric_limits<double>::quiet_NaN())
{
    if( tSF != 1.0 )
	// p(treatment seeking) is part of tree & handled generically, so we
	// don't have a way of modifying it.
	throw xml_scenario_error("treatment seeking heterogeneity not supported");
}
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
    // Note: we use Pathogenesis::COMPLICATED instead of Pathogenesis::SEVERE.
    Pathogenesis::State newState = pathogenesisModel->determineState (ageYears, withinHostModel);
    
    if ( Global::simulationTime == timeOfRecovery ) {
	if( pgState & Pathogenesis::DIRECT_DEATH ){
	    // Human dies this timestep (last day of risk of death)
	    _doomed = DOOMED_COMPLICATED;
	    
	    latestReport.update (Global::simulationTime, ageGroup, pgState);
	} else {
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
    }
    
    bool cmEvent = false;	// set true when we need to do case management
    if ( pgState & Pathogenesis::COMPLICATED ) {	// previous state: complicated
	// if severe, no events happen for course of medication
    } else if ( pgState & Pathogenesis::SICK ) {	// previous state: uncomplicated (UC/UC2/NMF)
	// if uc/uc2/nmf, the only event that can happen within 3 days is progression to severe/coinfection
	if ( newState & Pathogenesis::COMPLICATED ) {
	    pgState = Pathogenesis::State (pgState | newState);
	    cmEvent = true;
	}
    } else {	// previous state: healthy or delayed UC
	if ( newState & Pathogenesis::SICK ) {	// any (malarial/non-malarial) sickness
	    if( (pgState & Pathogenesis::PENDING_UC) == 0 ) {
		timeOfRecovery = Global::simulationTime + maxUCSeekingMemory;
		pgState = Pathogenesis::State (pgState | newState | Pathogenesis::PENDING_UC);
	    }
	}
	if( pgState & Pathogenesis::PENDING_UC ) {
	    if( random::uniform_01() < pImmediateUC)
		cmEvent = true;
	}
    }
    
    if ( cmEvent )	// note: new event can override pending delayed event
    {
	// If last treatment prescribed was in recent memory, consider second line.
	if (timeLastTreatment + Episode::healthSystemMemory > Global::simulationTime)
	    pgState = Pathogenesis::State (pgState | Pathogenesis::SECOND_CASE);
	
	caseStartTime = Global::simulationTime;
	
	if (pgState & Pathogenesis::MALARIA) {
	    if (util::ModelOptions::option (util::PENALISATION_EPISODES)) {
		withinHostModel.immunityPenalisation();
	    }
	}
	
	CMAuxOutput auxOut = ESCaseManagement::execute(
	    medicateQueue, pgState, withinHostModel, ageYears, ageGroup
	);
	
	if( medicateQueue.size() )	// I.E. some treatment was given (list is cleared either way)
	    timeLastTreatment = Global::simulationTime;
	
	if ( auxOut.hospitalisation != CMAuxOutput::NONE ) {	// in hospital
	    pgState = Pathogenesis::State (pgState | Pathogenesis::EVENT_IN_HOSPITAL);
	    
	    if( auxOut.hospitalisation == CMAuxOutput::DELAYED )
		caseStartTime++;
	}
	if ( auxOut.RDT_used ) {
	    Surveys.current->report_Clinical_RDTs (1);
	}
	
	// Case fatality rate (first day of illness)
	// P(death) is some fixed input scaled by age-specific CFR.
	if( (pgState & Pathogenesis::COMPLICATED)
	    && !(pgState & Pathogenesis::DIRECT_DEATH)
	) {
	    double pDeath = getPDeathInitial( ageYears );
	    // community fatality rate when not in hospital or delayed hospital entry
	    if( auxOut.hospitalisation != CMAuxOutput::IMMEDIATE )
		pDeath = (communityOddsMultiplier * pDeath) /
		    (1.0 + pDeath * (communityOddsMultiplier - 1.0) );
	    if (random::uniform_01() < pDeath) {
		pgState = Pathogenesis::State (pgState | Pathogenesis::DIRECT_DEATH | Pathogenesis::EVENT_FIRST_DAY);
		// Human is killed at end of time at risk
		timeOfRecovery += extraDaysAtRisk;
	    }
	    previousDensity = withinHostModel.getTotalDensity();
	}
    } else {
	// No new event (haven't changed state this timestep).
	
	// Case fatality rate (subsequent days)
	// Complicated case & at risk of death (note: extraDaysAtRisk <= 0)
	if( (pgState & Pathogenesis::COMPLICATED)
	    && !(pgState & Pathogenesis::DIRECT_DEATH)
	    && (Global::simulationTime < timeOfRecovery + extraDaysAtRisk)
	) {
	    // In complicated episodes, S(t), the probability of survival on
	    // subsequent days t, is described by log(S(t)) = -v(Y(t)/Y(t-1)),
	    // for parasite density Y(t). v_neg below is -v.
	    double parasiteReductionEffect = withinHostModel.getTotalDensity() / previousDensity;
	    double pDeath = 1.0 - exp( neg_v * parasiteReductionEffect );
	    // community fatality rate when not in hospital
	    if( !(pgState & Pathogenesis::EVENT_IN_HOSPITAL) )
		pDeath = (communityOddsMultiplier * pDeath) /
		(1.0 + pDeath * (communityOddsMultiplier - 1.0) );
	    if (random::uniform_01() < pDeath) {
		pgState = Pathogenesis::State (pgState | Pathogenesis::DIRECT_DEATH);
		// Human is killed at end of time at risk
		timeOfRecovery += extraDaysAtRisk;
	    }
	    previousDensity = withinHostModel.getTotalDensity();
	}
    }
    
    // Start of case. Not necessarily start of sickness due to treatment-seeking
    // delays and travel time.
    if( caseStartTime == Global::simulationTime ){
	// Patients in hospital are removed from the transmission cycle.
	// This should have an effect from the start of the next timestep.
	if (pgState & Pathogenesis::EVENT_IN_HOSPITAL)
	    hostTransmission.removeFromTransmission( true );
	
	if (pgState & Pathogenesis::COMPLICATED) {
	    // complicatedCaseDuration should to some respects be associated
	    // with medication duration, however ongoing medications after
	    // exiting hospital are OK and medications terminating before the
	    // end of hospitalization shouldn't matter too much if the person
	    // can't recieve new infections due to zero transmission in hospital.
	    timeOfRecovery = Global::simulationTime + complicatedCaseDuration;
	} else {
	    timeOfRecovery = Global::simulationTime + uncomplicatedCaseDuration;
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
    caseStartTime & stream;
    timeOfRecovery & stream;
    timeLastTreatment & stream;
    previousDensity & stream;
    medicateQueue & stream;
}
void ClinicalEventScheduler::checkpoint (ostream& stream) {
    ClinicalModel::checkpoint (stream);
    pgState & stream;
    caseStartTime & stream;
    timeOfRecovery & stream;
    timeLastTreatment & stream;
    previousDensity & stream;
    medicateQueue & stream;
}

} }