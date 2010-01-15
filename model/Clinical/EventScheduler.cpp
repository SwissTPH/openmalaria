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
	outcomes[it->getID()] = od;
    }
    outcomeMask = ev.getClinicalOutcomes().getSevereMask();
}


// -----  construction, destruction and checkpointing  -----

ClinicalEventScheduler::ClinicalEventScheduler (double cF, double tSF) :
        ClinicalModel (cF),
        pgState (Pathogenesis::NONE), pgChangeTimestep (Global::TIMESTEP_NEVER)
{}
ClinicalEventScheduler::~ClinicalEventScheduler() {}


// -----  other methods  -----

void ClinicalEventScheduler::doClinicalUpdate (WithinHost::WithinHostModel& withinHostModel, double ageYears)
{
    // Run pathogenesisModel
    // Note: we use Pathogenesis::COMPLICATED instead of Pathogenesis::SEVERE, though considering
    // the model's not designed to handle "coninfections" as presented by the Pathogenesis model
    // this is unimportant.
    Pathogenesis::State newState = pathogenesisModel->determineState (ageYears, withinHostModel);
    
    //TODO: changes to this condition (pending discussion)
    //new logic: if severe, nothing happens for "hospitalization days"
    // if uc/uc2/nmf, the only event that can happen within 3 days is severe
    
    /* Literally, if differences between (the combination of pgState and newState)
     * and pgState include SICK, MALARIA or COMPLICATED
     * or a second case of MALARIA and at least 3 days since last change */
    if ( ( ( (newState | pgState) ^ pgState) & (Pathogenesis::SICK | Pathogenesis::MALARIA | Pathogenesis::COMPLICATED)) ||
            (pgState & newState & Pathogenesis::MALARIA && Global::simulationTime >= pgChangeTimestep + 3))
    {
        if ( (newState & pgState) & Pathogenesis::MALARIA)
            newState = Pathogenesis::State (newState | Pathogenesis::SECOND_CASE);
        pgState = Pathogenesis::State (pgState | newState);
        SurveyAgeGroup ageGroup = ageYears;
        latestReport.update (Global::simulationTime, ageGroup, pgState);
        pgChangeTimestep = Global::simulationTime;
	
	if (pgState & Pathogenesis::MALARIA) {
	    if (util::ModelOptions::option (util::PENALISATION_EPISODES)) {
		withinHostModel.immunityPenalisation();
	    }
	}
	
        lastCmDecision = ESCaseManagement::execute (medicateQueue, pgState, withinHostModel, ageYears, ageGroup);
	
	// FIXME: set EVENT_IN_HOSPITAL from lastCmDecision
	// FIXME: report hospital entry
	//TODO report recoveries/death, and seq?
	
	if (pgState & Pathogenesis::COMPLICATED) {
	    //TODO: report sequelae?
	    // pgState = Pathogenesis::State (pgState | Pathogenesis::SEQUELAE);
	    
	    cmid id = lastCmDecision & outcomeMask;
	    OutcomeType::const_iterator oi = outcomes.find (id);
	    if (oi != outcomes.end ()) {
		if (gsl::rngUniform() < oi->second.pDeath) {
		    Surveys.current->reportHospitalizationDays (oi->second.hospitalizationDaysDeath);
		    
		    pgState = Pathogenesis::State (pgState | Pathogenesis::DIRECT_DEATH);
		    latestReport.update (Global::simulationTime, SurveyAgeGroup (ageYears), pgState);
		    //FIXME: kill after hospitalizationDaysDeath
		    _doomed = DOOMED_COMPLICATED; // kill human (removed from simulation next timestep)
		} else {
		    Surveys.current->reportHospitalizationDays (oi->second.hospitalizationDaysRecover);
		}
	    }
	    // FIXME: pgState still affects other things! Reset at some point?
	    //  pgState = Pathogenesis::NONE;
	}
    }


    if (pgState & Pathogenesis::INDIRECT_MORTALITY && _doomed == 0)
        _doomed = -Global::interval; // start indirect mortality countdown
    
    //FIXME: do we need to force time-outs now?
    if (pgState & Pathogenesis::SICK && latestReport.episodeEnd (Global::simulationTime)) {
        // Episode timeout: force recovery.
        // NOTE: An uncomplicated case occuring before this reset could be counted
        // UC2 but with treatment only occuring after this reset (when the case
        // should be counted UC) due to a treatment-seeking-delay. This can't be
        // corrected because the delay depends on the UC/UC2/etc. state leading to
        // a catch-22 situation, so DH, MP and VC decided to leave it like this.
        //TODO: also report EVENT_IN_HOSPITAL where relevant (patient _can_ be in a severe state)
        pgState = Pathogenesis::State (pgState | Pathogenesis::RECOVERY);
        latestReport.update (Global::simulationTime, SurveyAgeGroup (ageYears), pgState);
        pgState = Pathogenesis::NONE;
    }
    
    // Apply pending medications
    for (list<MedicateData>::iterator it = medicateQueue.begin(); it != medicateQueue.end();) {
        list<MedicateData>::iterator next = it;
        ++next;
        if (it->seekingDelay == 0) { // Medicate today's medications
            withinHostModel.medicate (it->abbrev, it->qty, it->time, ageYears);
            medicateQueue.erase (it);
            //TODO sort out reporting
        } else {   // and decrement treatment seeking delay for the rest
            it->seekingDelay--;
        }
        it = next;
    }
}


void ClinicalEventScheduler::checkpoint (istream& stream) {
    ClinicalModel::checkpoint (stream);
    int s;
    s & stream;
    pgState = Pathogenesis::State(s);
    pgChangeTimestep & stream;
    medicateQueue & stream;
    lastCmDecision & stream;
}
void ClinicalEventScheduler::checkpoint (ostream& stream) {
    ClinicalModel::checkpoint (stream);
    pgState & stream;
    pgChangeTimestep & stream;
    medicateQueue & stream;
    lastCmDecision & stream;
}

} }