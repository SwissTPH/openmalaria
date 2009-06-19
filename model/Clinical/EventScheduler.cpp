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
#include "GSLWrapper.h"
#include "WithinHostModel.h"
#include "simulation.h"

//TODO: move to XML
const int maxEpisodeLength = 28;

vector<double> ClinicalEventScheduler::caseManagementMaxAge;
vector<ClinicalEventScheduler::CaseManagementEndPoints> ClinicalEventScheduler::caseManagementEndPoints;


// -----  static init  -----

void ClinicalEventScheduler::init () {
  if (Global::interval != 1)
    throw xml_scenario_error("ClinicalEventScheduler is only designed for a 1-day timestep.");
  if (getCaseManagements() == NULL)
    throw xml_scenario_error ("ClinicalEventScheduler selected without caseManagements data in XML");
  
  const scnXml::CaseManagements::CaseManagementSequence& managements = getCaseManagements()->getCaseManagement();
  caseManagementMaxAge.resize (managements.size());
  caseManagementEndPoints.resize (managements.size());
  for (size_t i = 0; i < managements.size(); ++i) {
    caseManagementMaxAge[i] = managements[i].getMaxAgeYrs();
    CaseManagementEndPoints& endPoints = caseManagementEndPoints[i];
    endPoints.caseUC1 = readEndPoints (managements[i].getUc1());
    endPoints.caseUC2 = readEndPoints (managements[i].getUc2());
    endPoints.caseSev = readEndPoints (managements[i].getSev());
    endPoints.caseNMF = readEndPoints (managements[i].getNmf());
    
    const scnXml::Decisions::DecisionSequence& dSeq = managements[i].getDecisions().getDecision();
    for (scnXml::Decisions::DecisionSequence::const_iterator it = dSeq.begin(); it != dSeq.end(); ++it) {
      CaseTreatment ct;
      const scnXml::Decision::MedicateSequence& mSeq = it->getMedicate();
      ct.medications.resize (mSeq.size());
      for (size_t j = 0; j < mSeq.size(); ++j) {
	ct.medications[j].abbrev = mSeq[j].getName();
	ct.medications[j].qty = mSeq[j].getQty();
	ct.medications[j].delay = mSeq[j].getTime();
      }
      endPoints.decisions[it->getId()] = ct;
    }
  }
  if (caseManagementMaxAge[managements.size()-1] <= get_maximum_ageyrs())
    throw xml_scenario_error ("CaseManagements data doesn't cover all ages");
}

ClinicalEventScheduler::CaseTypeEndPoints ClinicalEventScheduler::readEndPoints (const scnXml::CaseType& caseType) {
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
    throw xml_scenario_error("EndPoint probabilities don't add up to 1.0 (CaseManagements)");
  // In any case, force it exactly 1.0 (because it could be slightly less,
  // meaning a random number x could have cumP<x<1.0, causing index errors.
  ret.cumProbs[caseTypeSeq.size()-1] = 1.0;
  return ret;
}


// -----  construction, destruction and checkpointing  -----

ClinicalEventScheduler::ClinicalEventScheduler (double cF, double tSF) :
    ClinicalModel (cF),
    pgState(Pathogenesis::NONE), reportState(Pathogenesis::NONE),
    pgChangeTimestep(TIMESTEP_NEVER), episodeStartTimestep(TIMESTEP_NEVER)
{}
ClinicalEventScheduler::~ClinicalEventScheduler() {
  // If something still to report, do so now
  if (Simulation::simulationTime > episodeStartTimestep + maxEpisodeLength) {
    Simulation::gMainSummary->report (reportState, _surveyPeriod, _ageGroup);
  }
}

ClinicalEventScheduler::ClinicalEventScheduler (istream& in) :
    ClinicalModel (in)
{
  int state;
  in >> state;
  pgState = (Pathogenesis::State) state;
  in >> state;
  reportState = (Pathogenesis::State) state;
  in >> pgChangeTimestep;
  in >> episodeStartTimestep;
  in >> _surveyPeriod;
  in >> _ageGroup;
}
void ClinicalEventScheduler::write (ostream& out) {
  pathogenesisModel->write (out);
  out << latestReport;
  out << _doomed << endl; 
  out << pgState << endl;
  out << reportState << endl;
  out << pgChangeTimestep << endl;
  out << episodeStartTimestep << endl;
  out << _surveyPeriod << endl;
  out << _ageGroup << endl;
}


// -----  other methods  -----

void ClinicalEventScheduler::doClinicalUpdate (WithinHostModel& withinHostModel, double ageYears) {
  // Run pathogenesisModel
  Pathogenesis::State newState = pathogenesisModel->determineState (ageYears, withinHostModel);
  
  /* Literally, if differences between (the combination of pgState and newState)
   * and pgState include SICK, MALARIA or COMPLICATED
   * or a second case of MALARIA and at least 3 days since last change */
  if ((((newState | pgState) ^ pgState) & (Pathogenesis::SICK | Pathogenesis::MALARIA | Pathogenesis::COMPLICATED)) ||
      (pgState & newState & Pathogenesis::MALARIA && Simulation::simulationTime >= pgChangeTimestep + 3))
  {
    // When an event occurs, if it's at least 28 days later than the first case,
    // we report the old episode and count the new case a new episode.
    if (Simulation::simulationTime > episodeStartTimestep + maxEpisodeLength) {
      Simulation::gMainSummary->report (reportState, _surveyPeriod, _ageGroup);
      episodeStartTimestep = Simulation::simulationTime;
      reportState = Pathogenesis::NONE;
      _surveyPeriod=Simulation::gMainSummary->getSurveyPeriod();
      _ageGroup=Simulation::gMainSummary->ageGroup(ageYears);
    }
    
    if ((newState & pgState) & Pathogenesis::MALARIA)
      newState = Pathogenesis::State (newState | Pathogenesis::SECOND_CASE);
    pgState = Pathogenesis::State (pgState | newState);
    reportState = Pathogenesis::State (reportState | pgState);
    pgChangeTimestep = Simulation::simulationTime;
    
    doCaseManagement (withinHostModel, ageYears);
  }
  
  
  if (pgState & Pathogenesis::INDIRECT_MORTALITY && _doomed == 0)
    _doomed = -Global::interval;	// start indirect mortality countdown
  
  
  if (pgState & Pathogenesis::COMPLICATED) {
    if (Simulation::simulationTime >= pgChangeTimestep + 10) {
      // force recovery after 10 days
      if (W_UNIFORM() < 0.02)
	reportState = Pathogenesis::State (reportState | Pathogenesis::SEQUELAE);
      pgState = Pathogenesis::NONE;
    } else {
      //TODO: insert correct probabilities
      const double pRecover = 0.1;
      const double pSequelae = 0.02;
      const double pDeath = 0.03;
      double rand = W_UNIFORM();
      if (rand < pRecover) {
	if (rand < pSequelae*pRecover && Simulation::simulationTime >= pgChangeTimestep + 5)
	  reportState = Pathogenesis::State (reportState | Pathogenesis::SEQUELAE);
	pgState = Pathogenesis::NONE;
      } else if (rand < pRecover+pDeath) {
	reportState = Pathogenesis::State (reportState | Pathogenesis::DIRECT_DEATH);
	_doomed = 4;	// kill human (removed from simulation next timestep)
      }
      // else stay in this state
    }
  } else if (Simulation::simulationTime >= episodeStartTimestep + maxEpisodeLength) {
    // End of what's counted as episode. We only do reporting on death or the
    // next event.
    // NOTE: An uncomplicated case occuring before this reset could be counted
    // UC2 but with treatment only occuring after this reset (when the case
    // should be counted UC) due to a treatment-seeking-delay. This can't be
    // corrected because the delay depends on the UC/UC2/etc. state leading to
    // a catch-22 situation, so DH, MP and VC decided to leave it like this.
    pgState = Pathogenesis::NONE;
  }
  
  for (list<MedicateData>::iterator it = medicateQueue.begin(); it != medicateQueue.end();) {
    list<MedicateData>::iterator next = it;
    ++next;
    if (it->delay < 24) {	// Medicate today's medications
      withinHostModel.medicate(it->abbrev, it->qty, it->delay);
      medicateQueue.erase(it);
      //TODO sort out reporting
    } else {			// and decrement delay for the rest
      it->delay -= 24;
    }
    it = next;
  }
}

void ClinicalEventScheduler::doCaseManagement (WithinHostModel& withinHostModel, double ageYears) {
#ifdef DEBUG
  if (!(pgState & Pathogenesis::SICK))
    throw domain_error("doCaseManagement shouldn't be called if not sick");
#endif
  
  // We always remove any queued medications.
  medicateQueue.clear();
  
  size_t ageIndex = 0;
  while (ageYears > caseManagementMaxAge[ageIndex]) {
    ++ageIndex;
#ifdef DEBUG
    if (ageIndex >= caseManagementMaxAge.size()) {
      ostringstream x;
      x << "Individual's age (" << ageYears << ") is over maximum age which has caseManagement data in XML (" << CaseManagementMaxAge[ageIndex].first << ")";
      throw xml_scenario_error(x.str());
    }
#endif
  }
  
  const CaseTypeEndPoints* endPoints;
  if (pgState & Pathogenesis::MALARIA) {
    if (pgState & Pathogenesis::COMPLICATED)
      endPoints = &caseManagementEndPoints[ageIndex].caseSev;
    else if (pgState & Pathogenesis::SECOND_CASE)
      endPoints = &caseManagementEndPoints[ageIndex].caseUC2;
    else
      endPoints = &caseManagementEndPoints[ageIndex].caseUC1;
  } else /*if (pgState & Pathogenesis::SICK) [true by above check]*/ {	// sick but not from malaria
    endPoints = &caseManagementEndPoints[ageIndex].caseNMF;
  }
  
  double randCum = W_UNIFORM();
  size_t decisionIndex = 0;
  while (endPoints->cumProbs[decisionIndex] < randCum)
    ++decisionIndex;
  
  CaseTreatment& decision = caseManagementEndPoints[ageIndex].decisions[endPoints->decisions[decisionIndex]];
  for (vector<MedicateData>::iterator it = decision.medications.begin(); it != decision.medications.end(); ++it)
    medicateQueue.push_back (*it);
}
