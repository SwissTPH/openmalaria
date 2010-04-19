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
#include "Surveys.h"
#include "util/ModelOptions.hpp"
#include "util/errors.hpp"

namespace OM { namespace Clinical {
    using namespace OM::util;

int ClinicalEventScheduler::uncomplicatedCaseDuration;
int ClinicalEventScheduler::complicatedCaseDuration;
double ClinicalEventScheduler::pDeathInitial;
double ClinicalEventScheduler::parasiteReductionEffectiveness;

//FIXME: age correction factor
double pDeathAgeFactor (double ageYears) {
    return 1.0;
}


// -----  static init  -----

void ClinicalEventScheduler::init ()
{
    if (Global::interval != 1)
        throw util::xml_scenario_error ("ClinicalEventScheduler is only designed for a 1-day timestep.");
    if (! (util::ModelOptions::option (util::INCLUDES_PK_PD)))
        throw util::xml_scenario_error ("ClinicalEventScheduler requires INCLUDES_PK_PD");
    
    // We assume EventScheduler data was already confirmed present:
    const scnXml::HSEventScheduler& ev = InputData().getHealthSystem().getEventScheduler().get();
    
    ESCaseManagement::init ();
    /*
    const scnXml::ClinicalOutcomes::OutcomesSequence& ocSeq = ev.getClinicalOutcomes().getOutcomes();
    for (scnXml::ClinicalOutcomes::OutcomesConstIterator it = ocSeq.begin(); it != ocSeq.end(); ++it) {
	OutcomeData od;
	od.pDeath = it->getPDeath();
	od.hospitalizationDaysDeath = it->getHospitalizationDaysDeath ();
	od.hospitalizationDaysRecover = it->getHospitalizationDaysRecover ();
	//FIXME: assert hospitalizationDaysRecover equals length of medication course in days
	outcomes[it->getID()] = od;
    }
    outcomeMask = ev.getClinicalOutcomes().getSevereMask();*/
    
    //FIXME: initialize properly
    uncomplicatedCaseDuration = 3;
    complicatedCaseDuration = 6;
    pDeathInitial = 0.2;
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

void ClinicalEventScheduler::doClinicalUpdate (WithinHost::WithinHostModel& withinHostModel, double ageYears, SurveyAgeGroup ageGroup)
{
    // Run pathogenesisModel
    // Note: we use Pathogenesis::COMPLICATED instead of Pathogenesis::SEVERE, though considering
    // the model's not designed to handle "coninfections" as presented by the Pathogenesis model
    // this is unimportant.
    Pathogenesis::State newState = pathogenesisModel->determineState (ageYears, withinHostModel);
    
    if ( Global::simulationTime == timeOfRecovery ) {
	/*{//FIXME: report recovery/sequelae
	    // Note: setting SEQUELAE / RECOVERY here only affects reporting, not model
	    //TODO: report sequelae?
	    // if ( Sequelae )
	    //     pgState = Pathogenesis::State (pgState | Pathogenesis::SEQUELAE);
	    // else
	    pgState = Pathogenesis::State (pgState | Pathogenesis::RECOVERY);
	}*/
	
	// Individual recovers (and is immediately susceptible to new cases)
	pgState = Pathogenesis::NONE;	// recovery (reset to healthy state)
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
        latestReport.update (Global::simulationTime, ageGroup, pgState);
	
	if (pgState & Pathogenesis::MALARIA) {
	    if (util::ModelOptions::option (util::PENALISATION_EPISODES)) {
		withinHostModel.immunityPenalisation();
	    }
	}
	
	//TODO: get quality of case management and work into p(death)
	ESCaseManagement::execute (medicateQueue, pgState, withinHostModel, ageYears, ageGroup);
	
	if (false /*FIXME: RDT used*/) {
	    Surveys.current->report_Clinical_RDTs (1);
	}
	
	if (pgState & Pathogenesis::COMPLICATED) {
	    //TODO: zero infectiousness when in hospital
	    
	    //TODO: document (first day of case fatality model)
	    if (random::uniform_01() < pDeathInitial) {
		pgState = Pathogenesis::State (pgState | Pathogenesis::DIRECT_DEATH);
		latestReport.update (Global::simulationTime, ageGroup, pgState);
		
		// Human dies this timestep; we still allow medication of
		// today's treatments for costing.
		_doomed = DOOMED_COMPLICATED;
	    }
	    previousDensity = withinHostModel.getTotalDensity();
	    
	    // Check: is patient in hospital? Reporting only.
// 	    if ( lastCmDecision & (Decision::TREATMENT_HOSPITAL | Decision::TREATMENT_DEL_HOSPITAL) ) {	// in hospital
// 		pgState = Pathogenesis::State (pgState | Pathogenesis::EVENT_IN_HOSPITAL);
// 		Surveys.current->report_Clinical_HospitalEntries (1);
// 		Surveys.current->report_Clinical_HospitalizationDays (medicationDuration);
// 	    }
	    
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
	
	//TODO: document (subsequent days of case fatality model)
	if (pgState & Pathogenesis::COMPLICATED) {
	    double parasiteReductionEffect = pow(
		previousDensity / (previousDensity - withinHostModel.getTotalDensity()),
		parasiteReductionEffectiveness
	    );
	    double pDeath = 1.0 - exp(-pDeathAgeFactor(ageYears) * parasiteReductionEffect);
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
    uncomplicatedCaseDuration & stream;
    complicatedCaseDuration & stream;
    pDeathInitial & stream;
    parasiteReductionEffectiveness & stream;
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
    uncomplicatedCaseDuration & stream;
    complicatedCaseDuration & stream;
    pDeathInitial & stream;
    parasiteReductionEffectiveness & stream;
    pgState & stream;
    lastEventTime & stream;
    timeOfRecovery & stream;
    previousDensity & stream;
    medicateQueue & stream;
}

} }