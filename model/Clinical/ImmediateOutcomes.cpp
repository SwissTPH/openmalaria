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

#include "Clinical/ImmediateOutcomes.h"


namespace OM { namespace Clinical {

// -----  static init  -----

void ClinicalImmediateOutcomes::initParameters () {
  OldCaseManagement::init();
}


// -----  construction and destruction  -----

ClinicalImmediateOutcomes::ClinicalImmediateOutcomes (double cF, double tSF) :
    ClinicalModel (cF),
    caseManagement(new OldCaseManagement (tSF))
{}
ClinicalImmediateOutcomes::~ClinicalImmediateOutcomes() {
  delete caseManagement; 
}


// -----  other methods  -----

void ClinicalImmediateOutcomes::massDrugAdministration(OM::WithinHost::WithinHostModel& withinHostModel, double ageYears) {
    /* TODO: here we assume a 100% clearance rate for the MDA drug we use.
       This is not consistent with the way we treat according to the Health
       system description. The default clearance rate for MDA should be 100%
       since this simulates what was meant to happen in Garki.  We can change
       this by introducing an optional clearance rate that can be < 100%.
    */
    // We need to pass the is-severe state for the IPT code.
    withinHostModel.clearInfections(latestReport.getState() == Pathogenesis::STATE_SEVERE);
}

void ClinicalImmediateOutcomes::doClinicalUpdate (WithinHost::WithinHostModel& withinHostModel, double ageYears, SurveyAgeGroup ageGroup) {
  caseManagement->doCaseManagement (pathogenesisModel->determineState (ageYears, withinHostModel),
				    withinHostModel,
				    latestReport,
				    ageYears,
				    ageGroup,
				    _doomed);
}


void ClinicalImmediateOutcomes::checkpoint (istream& stream) {
    ClinicalModel::checkpoint (stream);
    (*caseManagement) & stream;
}
void ClinicalImmediateOutcomes::checkpoint (ostream& stream) {
    ClinicalModel::checkpoint (stream);
    (*caseManagement) & stream;
}

} }