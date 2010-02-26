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
#include "util/gsl.h"
#include "WithinHost/WithinHostModel.h"
#include "Surveys.h"
#include "util/ModelOptions.hpp"
#include "util/errors.hpp"

namespace OM { namespace Clinical {

ClinicalEventScheduler::OutcomeType ClinicalEventScheduler::outcomes;
cmid ClinicalEventScheduler::outcomeMask;


// -----  static init  -----

void ClinicalEventScheduler::init ()
{
    if (Global::interval != 1)
        throw util::xml_scenario_error ("ClinicalEventScheduler is only designed for a 1-day timestep.");
    if (! (util::ModelOptions::option (util::INCLUDES_PK_PD)))
        throw util::xml_scenario_error ("ClinicalEventScheduler requires INCLUDES_PK_PD");
    
    const scnXml::EventScheduler& ev = InputData.getEventScheduler();
    Episode::healthSystemMemory = ev.getHealthSystemMemory();
    
    ESCaseManagement::init ();
    
    const scnXml::ClinicalOutcomes::OutcomesSequence& ocSeq = ev.getClinicalOutcomes().getOutcomes();
    for (scnXml::ClinicalOutcomes::OutcomesConstIterator it = ocSeq.begin(); it != ocSeq.end(); ++it) {
	OutcomeData od;
	od.pDeath = it->getPDeath();
	od.hospitalizationDaysDeath = it->getHospitalizationDaysDeath ();
	od.hospitalizationDaysRecover = it->getHospitalizationDaysRecover ();
	//FIXME: assert hospitalizationDaysRecover equals length of medication course in days
	outcomes[it->getID()] = od;
    }
    outcomeMask = ev.getClinicalOutcomes().getSevereMask();
}


// -----  construction, destruction and checkpointing  -----

ClinicalEventScheduler::ClinicalEventScheduler (double cF, double tSF) :
        ClinicalModel (cF),
        pgState (Pathogenesis::NONE), timeHealthyOrDead (Global::TIMESTEP_NEVER)
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

void ClinicalEventScheduler::doClinicalUpdate (WithinHost::WithinHostModel& withinHostModel, double ageYears)
{
    // Run pathogenesisModel
    // Note: we use Pathogenesis::COMPLICATED instead of Pathogenesis::SEVERE, though considering
    // the model's not designed to handle "coninfections" as presented by the Pathogenesis model
    // this is unimportant.
    Pathogenesis::State newState = pathogenesisModel->determineState (ageYears, withinHostModel);
    
    if ( Global::simulationTime == timeHealthyOrDead ) {
	if ( pgState & Pathogenesis::DIRECT_DEATH ) {
	    // Note: killed now, no further clinical events but still a source of infection to mosquitoes until next timestep
	    _doomed = DOOMED_COMPLICATED; // kill human (removed from simulation next timestep)
	    return;
	}
	
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
        SurveyAgeGroup ageGroup = ageYears;
        latestReport.update (Global::simulationTime, ageGroup, pgState);
	
	if (pgState & Pathogenesis::MALARIA) {
	    if (util::ModelOptions::option (util::PENALISATION_EPISODES)) {
		withinHostModel.immunityPenalisation();
	    }
	}
	
	// Decision ID of last case management run
	cmid lastCmDecision = ESCaseManagement::execute (medicateQueue, pgState, withinHostModel, ageYears, ageGroup);
	
	if ( lastCmDecision & Decision::TEST_RDT ) {
	    Surveys.current->report_Clinical_RDTs (1);
	}
	
	if (pgState & Pathogenesis::COMPLICATED) {
	    // Find outcome probabilities/durations associated with CM-tree path
	    cmid id = lastCmDecision & outcomeMask;
	    OutcomeType::const_iterator oi = outcomes.find (id);
	    if (oi == outcomes.end ()) {
		ostringstream msg;
		msg << "No 'Outcome' data for cmid " << id;
		throw util::xml_scenario_error(msg.str());
	    }
	    
	    int medicationDuration;
	    if (rng::uniform01() < oi->second.pDeath) {
		medicationDuration = oi->second.hospitalizationDaysDeath;
		
		pgState = Pathogenesis::State (pgState | Pathogenesis::DIRECT_DEATH);
	    } else {
		medicationDuration = oi->second.hospitalizationDaysRecover;
		
		// Note: setting SEQUELAE / RECOVERY here only affects reporting, not model
		//TODO: report sequelae?
		// if ( Sequelae )
		//     pgState = Pathogenesis::State (pgState | Pathogenesis::SEQUELAE);
		// else
		pgState = Pathogenesis::State (pgState | Pathogenesis::RECOVERY);
	    }
	    
	    // Check: is patient in hospital? Reporting only.
	    if ( lastCmDecision & (Decision::TREATMENT_HOSPITAL | Decision::TREATMENT_DEL_HOSPITAL) ) {	// in hospital
		pgState = Pathogenesis::State (pgState | Pathogenesis::EVENT_IN_HOSPITAL);
		Surveys.current->report_Clinical_HospitalEntries (1);
		Surveys.current->report_Clinical_HospitalizationDays (medicationDuration);
	    }
	    
	    // Report: recovery/seq./death, in/out of hospital
	    latestReport.update (Global::simulationTime, SurveyAgeGroup (ageYears), pgState);
	    timeHealthyOrDead = Global::simulationTime + medicationDuration;
	} else {
	    timeHealthyOrDead = Global::simulationTime + 3;
	}
    }


    if (pgState & Pathogenesis::INDIRECT_MORTALITY && _doomed == 0)
        _doomed = -Global::interval; // start indirect mortality countdown
    
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
    timeHealthyOrDead & stream;
    medicateQueue & stream;
}
void ClinicalEventScheduler::checkpoint (ostream& stream) {
    ClinicalModel::checkpoint (stream);
    pgState & stream;
    timeHealthyOrDead & stream;
    medicateQueue & stream;
}

} }