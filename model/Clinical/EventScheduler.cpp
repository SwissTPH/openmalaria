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

double ClinicalEventScheduler::pDeathTable[TREATMENT_NUM_TYPES * PTABLE_NUM_DAYS];
double ClinicalEventScheduler::pRecoverTable[TREATMENT_NUM_TYPES * PTABLE_NUM_DAYS];
vector<double> ClinicalEventScheduler::caseManagementMaxAge;
vector<ClinicalEventScheduler::CaseManagementEndPoints> ClinicalEventScheduler::caseManagementEndPoints;


// -----  static init  -----

void ClinicalEventScheduler::init ()
{
  if (Global::interval != 1)
    throw xml_scenario_error ("ClinicalEventScheduler is only designed for a 1-day timestep.");
  if (! (Global::modelVersion & INCLUDES_PK_PD))
    throw xml_scenario_error ("ClinicalEventScheduler requires INCLUDES_PK_PD");
  if (getCaseManagements() == NULL)
    throw xml_scenario_error ("ClinicalEventScheduler selected without caseManagements data in XML");

  Episode::reportingPeriodMemory = getCaseManagements()->getReportingPeriodMemory();
  
  //TODO: set properly from XML; also probability of sequelae
  for (size_t t = 0; t < TREATMENT_NUM_TYPES * PTABLE_NUM_DAYS; ++t) {
    pDeathTable[t] = 0.03;
    pRecoverTable[t] = 0.1;
  }
  
  const scnXml::CaseManagements::CaseManagementSequence& managements = getCaseManagements()->getCaseManagement();
  caseManagementMaxAge.resize (managements.size());
  caseManagementEndPoints.resize (managements.size());
  for (size_t i = 0; i < managements.size(); ++i) {
    caseManagementMaxAge[i] = managements[i].getMaxAgeYrs();
    CaseManagementEndPoints& endPoints = caseManagementEndPoints[i];
    endPoints.caseUC1 = readEndPoints (managements[i].getUc1());
    endPoints.caseUC2 = readEndPoints (managements[i].getUc2());
    endPoints.caseSev = readEndPoints (managements[i].getSev());
    endPoints.caseNMFWithParasites = readEndPoints (managements[i].getNmfP());
    endPoints.caseNMFWithoutParasites = readEndPoints (managements[i].getNmfNP());

    const scnXml::Decisions::DecisionSequence& dSeq = managements[i].getDecisions().getDecision();
    for (scnXml::Decisions::DecisionSequence::const_iterator it = dSeq.begin(); it != dSeq.end(); ++it) {
      CaseTreatment ct;
      const scnXml::Decision::MedicateSequence& mSeq = it->getMedicate();
      ct.medications.resize (mSeq.size());
      for (size_t j = 0; j < mSeq.size(); ++j) {
        ct.medications[j].abbrev = mSeq[j].getName();
        ct.medications[j].qty = mSeq[j].getQty();
        ct.medications[j].time = mSeq[j].getTime();
      }
      endPoints.decisions[it->getId() ] = ct;
    }
  }
  if (caseManagementMaxAge[managements.size()-1] <= get_maximum_ageyrs())
    throw xml_scenario_error ("CaseManagements data doesn't cover all ages");
}

ClinicalEventScheduler::CaseTypeEndPoints ClinicalEventScheduler::readEndPoints (const scnXml::CaseType& caseType)
{
  CaseTypeEndPoints ret;
  const scnXml::CaseType::EndPointSequence& caseTypeSeq = caseType.getEndPoint();
  ret.cumProbs.resize (caseTypeSeq.size());
  ret.decisions.resize (caseTypeSeq.size());
  double cumP = 0.0;
  for (size_t i = 0; i < caseTypeSeq.size(); ++i) {
    cumP += caseTypeSeq[i].getP();
    ret.cumProbs[i] = cumP;
    ret.decisions[i] = caseTypeSeq[i].getDecision();
  }
  // Test cumP is approx. 1.0 (in case the XML is wrong).
  if (cumP < 0.999 || cumP > 1.001)
    throw xml_scenario_error ("EndPoint probabilities don't add up to 1.0 (CaseManagements)");
  // In any case, force it exactly 1.0 (because it could be slightly less,
  // meaning a random number x could have cumP<x<1.0, causing index errors.
  ret.cumProbs[caseTypeSeq.size()-1] = 1.0;
  return ret;
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
    latestReport.update (Simulation::simulationTime, Simulation::gMainSummary->ageGroup (ageYears), pgState);
    pgChangeTimestep = Simulation::simulationTime;

    doCaseManagement (withinHostModel, ageYears);
  }


  if (pgState & Pathogenesis::INDIRECT_MORTALITY && _doomed == 0)
    _doomed = -Global::interval; // start indirect mortality countdown

  //TODO: Also call immunityPenalisation() like previous model?

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
      latestReport.update (Simulation::simulationTime, Simulation::gMainSummary->ageGroup (ageYears), pgState);
      pgState = Pathogenesis::NONE;
    } else if (daySinceSevere >= 1) {	//TODO: do we delay one day?
      // determine case fatality rates for day1, day2, day3 (remaining days are at day 3 probabilities)
      int delayIndex = daySinceSevere-1;	//TODO: is this right?
      if (delayIndex > 2) delayIndex = 2;	// use 3rd-day's value for any later days
      int index = (lastCmDecision & TREATMENT_MASK) >> TREATMENT_SHIFT;
      if (index >= TREATMENT_NUM_TYPES) throw runtime_error ("CM returned invalid treatment type");
      index += delayIndex * TREATMENT_NUM_TYPES;
      double pDeath = pDeathTable[index];
      double pRecover = pRecoverTable[index];
      
      double rand = gsl::rngUniform();
      if (rand < pRecover) {
        if (rand < pSequelae*pRecover && Simulation::simulationTime >= pgChangeTimestep + 5)
          pgState = Pathogenesis::State (pgState | Pathogenesis::SEQUELAE);
        else
          pgState = Pathogenesis::State (pgState | Pathogenesis::RECOVERY);
        latestReport.update (Simulation::simulationTime, Simulation::gMainSummary->ageGroup (ageYears), pgState);
        pgState = Pathogenesis::NONE;
      } else if (rand < pRecover + pDeath) {
        pgState = Pathogenesis::State (pgState | Pathogenesis::DIRECT_DEATH);
        latestReport.update (Simulation::simulationTime, Simulation::gMainSummary->ageGroup (ageYears), pgState);
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
    latestReport.update (Simulation::simulationTime, Simulation::gMainSummary->ageGroup (ageYears), pgState);
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

void ClinicalEventScheduler::doCaseManagement (WithinHostModel& withinHostModel, double ageYears)
{
#ifndef NDEBUG
  if (! (pgState & Pathogenesis::SICK))
    throw domain_error ("doCaseManagement shouldn't be called if not sick");
#endif

  // We always remove any queued medications.
  medicateQueue.clear();

  size_t ageIndex = 0;
  while (ageYears > caseManagementMaxAge[ageIndex]) {
    ++ageIndex;
#ifndef NDEBUG
    if (ageIndex >= caseManagementMaxAge.size()) {
      ostringstream x;
      x << "Individual's age (" << ageYears << ") is over maximum age which has caseManagement data in XML (" << caseManagementMaxAge[ageIndex-1] << ")";
      throw xml_scenario_error (x.str());
    }
#endif
  }

  const CaseTypeEndPoints* endPoints;
  if (pgState & Pathogenesis::MALARIA) { // NOTE: report treatment shouldn't be done like this so it's handled correctly when treatment is cancelled
    if (pgState & Pathogenesis::COMPLICATED) {
      endPoints = &caseManagementEndPoints[ageIndex].caseSev;
      Simulation::gMainSummary->reportTreatment (Simulation::gMainSummary->ageGroup (ageYears), 3);
    } else if (pgState & Pathogenesis::SECOND_CASE) {
      endPoints = &caseManagementEndPoints[ageIndex].caseUC2;
      Simulation::gMainSummary->reportTreatment (Simulation::gMainSummary->ageGroup (ageYears), 2);
    } else {
      endPoints = &caseManagementEndPoints[ageIndex].caseUC1;
      Simulation::gMainSummary->reportTreatment (Simulation::gMainSummary->ageGroup (ageYears), 1);
    }
  } else /*if (pgState & Pathogenesis::SICK) [true by above check]*/ { // sick but not from malaria
    if (withinHostModel.parasiteDensityDetectible())
      endPoints = &caseManagementEndPoints[ageIndex].caseNMFWithParasites;
    else
      endPoints = &caseManagementEndPoints[ageIndex].caseNMFWithoutParasites;
  }

  double randCum = gsl::rngUniform();
  size_t decisionIndex = 0;
  while (endPoints->cumProbs[decisionIndex] < randCum)
    ++decisionIndex;

  CaseTreatment& decision = caseManagementEndPoints[ageIndex].decisions[endPoints->decisions[decisionIndex]];
  for (vector<MedicateData>::iterator it = decision.medications.begin(); it != decision.medications.end(); ++it) {
    medicateQueue.push_back (*it);
    medicateQueue.back().seekingDelay = endPoints->decisions[decisionIndex] % 10; // last digit
    lastCmDecision = endPoints->decisions[decisionIndex];
  }
}
