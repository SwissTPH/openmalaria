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
#include "Simulation.h"
#include "Surveys.h"


namespace Clinical {


// -----  static init  -----

void ClinicalEventScheduler::init ()
{
    if (Global::interval != 1)
        throw xml_scenario_error ("ClinicalEventScheduler is only designed for a 1-day timestep.");
    if (! (Global::modelVersion & INCLUDES_PK_PD))
        throw xml_scenario_error ("ClinicalEventScheduler requires INCLUDES_PK_PD");

    Episode::healthSystemMemory = getEventScheduler().getHealthSystemMemory();

    //TODO: set properly from XML; also probability of sequelae
    /*for (size_t t = 0; t < PTABLE_SIZE; ++t) {
        pDeathTable[t] = 0.03;
        pRecoverTable[t] = 0.1;
    }*/

    ESCaseManagement::init ();
}


// -----  construction, destruction and checkpointing  -----

ClinicalEventScheduler::ClinicalEventScheduler (double cF, double tSF) :
        ClinicalModel (cF),
        pgState (Pathogenesis::NONE), pgChangeTimestep (TIMESTEP_NEVER)
{}
ClinicalEventScheduler::~ClinicalEventScheduler() {}

ClinicalEventScheduler::ClinicalEventScheduler (istream& in) :
        ClinicalModel (in)
{
    int x;
    in >> x;
    pgState = (Pathogenesis::State) x;
    in >> pgChangeTimestep;
    in >> x;
    for (; x > 0; --x) {
        MedicateData md;
        in >> md.abbrev;
        in >> md.qty;
        in >> md.time;
        in >> md.seekingDelay;
        medicateQueue.push_back (md);
    }
    in >> lastCmDecision;
}
void ClinicalEventScheduler::write (ostream& out)
{
    ClinicalModel::write (out);
    out << pgState << endl;
    out << pgChangeTimestep << endl;
    out << medicateQueue.size() << endl;
    for (list<MedicateData>::iterator i = medicateQueue.begin(); i != medicateQueue.end(); ++i) {
        out << i->abbrev << endl;
        out << i->qty << endl;
        out << i->time << endl;
        out << i->seekingDelay << endl;
    }
    out << lastCmDecision << endl;
}


// -----  other methods  -----

void ClinicalEventScheduler::doClinicalUpdate (WithinHostModel& withinHostModel, double ageYears)
{
    // Run pathogenesisModel
    Pathogenesis::State newState = pathogenesisModel->determineState (ageYears, withinHostModel);

    /* Literally, if differences between (the combination of pgState and newState)
     * and pgState include SICK, MALARIA or COMPLICATED
     * or a second case of MALARIA and at least 3 days since last change */
    if ( ( ( (newState | pgState) ^ pgState) & (Pathogenesis::SICK | Pathogenesis::MALARIA | Pathogenesis::COMPLICATED)) ||
            (pgState & newState & Pathogenesis::MALARIA && Simulation::simulationTime >= pgChangeTimestep + 3)) {
        if ( (newState & pgState) & Pathogenesis::MALARIA)
            newState = Pathogenesis::State (newState | Pathogenesis::SECOND_CASE);
        pgState = Pathogenesis::State (pgState | newState);
        SurveyAgeGroup ageGroup = ageYears;
        latestReport.update (Simulation::simulationTime, ageGroup, pgState);
        pgChangeTimestep = Simulation::simulationTime;

        ESCaseManagement::execute (medicateQueue, pgState, withinHostModel, ageYears, ageGroup);
    }


    if (pgState & Pathogenesis::INDIRECT_MORTALITY && _doomed == 0)
        _doomed = -Global::interval; // start indirect mortality countdown

    if (pgState & Pathogenesis::MALARIA) {
        if (Global::modelVersion & PENALISATION_EPISODES) {
            withinHostModel.immunityPenalisation();
        }
    }

    if (pgState & Pathogenesis::COMPLICATED) {
        //TODO: Also set Pathogenesis::EVENT_IN_HOSPITAL where relevant:
        //reportState = Pathogenesis::State (reportState | Pathogenesis::EVENT_IN_HOSPITAL);
        const double pSequelae = 0.02; //prob of sequelae is constant of recovereds
        // TODO pSequelae should come from xml

        int daySinceSevere = Simulation::simulationTime - pgChangeTimestep;

        if (daySinceSevere >= 10) {
            // force recovery after 10 days
            if (gsl::rngUniform() < pSequelae)
                pgState = Pathogenesis::State (pgState | Pathogenesis::SEQUELAE);
            else
                pgState = Pathogenesis::State (pgState | Pathogenesis::RECOVERY);
            latestReport.update (Simulation::simulationTime, SurveyAgeGroup (ageYears), pgState);
            pgState = Pathogenesis::NONE;
        } else if (daySinceSevere >= 1) { //TODO: do we delay one day?
            // determine case fatality rates for day1, day2, day3 (remaining days are at day 3 probabilities)
            int delayIndex = daySinceSevere - 1; //TODO: is this right?
            if (delayIndex > 2) delayIndex = 2; // use 3rd-day's value for any later days
            // NOTE: delayIndex removed from lookup since it's not in VC's tables

            /*FIXME
	    int index = (lastCmDecision & TREATMENT_MASK) >> TREATMENT_SHIFT;
            int drugIndex = lastCmDecision & DRUG_MASK;
            if (index >= TREATMENT_NUM_TYPES || drugIndex < DRUG_FIRST_SEV) {
                stringstream msg;
                msg << "CM decision invalid: 0x" << hex << lastCmDecision;
                throw runtime_error (msg.str());
            }
            drugIndex = (drugIndex - DRUG_FIRST_SEV) >> DRUG_SHIFT;
            index += drugIndex * TREATMENT_NUM_TYPES;*/
	    
            double pDeath = 0.03;	//pDeathTable[index];
            double pRecover = 0.1;	//pRecoverTable[index];

            double rand = gsl::rngUniform();
            if (rand < pRecover) {
                if (rand < pSequelae*pRecover && Simulation::simulationTime >= pgChangeTimestep + 5)
                    pgState = Pathogenesis::State (pgState | Pathogenesis::SEQUELAE);
                else
                    pgState = Pathogenesis::State (pgState | Pathogenesis::RECOVERY);
                latestReport.update (Simulation::simulationTime, SurveyAgeGroup (ageYears), pgState);
                pgState = Pathogenesis::NONE;
            } else if (rand < pRecover + pDeath) {
                pgState = Pathogenesis::State (pgState | Pathogenesis::DIRECT_DEATH);
                latestReport.update (Simulation::simulationTime, SurveyAgeGroup (ageYears), pgState);
                _doomed = DOOMED_COMPLICATED; // kill human (removed from simulation next timestep)
            }
            // else stay in this state
        }
    } else if (pgState & Pathogenesis::SICK && latestReport.episodeEnd (Simulation::simulationTime)) {
        // Episode timeout: force recovery.
        // NOTE: An uncomplicated case occuring before this reset could be counted
        // UC2 but with treatment only occuring after this reset (when the case
        // should be counted UC) due to a treatment-seeking-delay. This can't be
        // corrected because the delay depends on the UC/UC2/etc. state leading to
        // a catch-22 situation, so DH, MP and VC decided to leave it like this.
        //TODO: also report EVENT_IN_HOSPITAL where relevant (patient _can_ be in a severe state)
        pgState = Pathogenesis::State (pgState | Pathogenesis::RECOVERY);
        latestReport.update (Simulation::simulationTime, SurveyAgeGroup (ageYears), pgState);
        pgState = Pathogenesis::NONE;
    }

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

}