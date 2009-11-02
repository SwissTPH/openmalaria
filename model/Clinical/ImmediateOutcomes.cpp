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


// -----  static init  -----

void ClinicalImmediateOutcomes::initParameters () {
  Episode::reportingPeriodMemory = OldCaseManagement::init();
}


// -----  construction, destruction and checkpointing  -----

ClinicalImmediateOutcomes::ClinicalImmediateOutcomes (double cF, double tSF) :
    ClinicalModel (cF),
    caseManagement(new OldCaseManagement (tSF))
{}
ClinicalImmediateOutcomes::~ClinicalImmediateOutcomes() {
  delete caseManagement; 
}

ClinicalImmediateOutcomes::ClinicalImmediateOutcomes (istream& in) :
    ClinicalModel (in),
    caseManagement(new OldCaseManagement (in))
{}
void ClinicalImmediateOutcomes::write (ostream& out) {
  ClinicalModel::write (out);
  caseManagement->write (out);
}


// -----  other methods  -----

void ClinicalImmediateOutcomes::doClinicalUpdate (WithinHostModel& withinHostModel, double ageYears) {
  caseManagement->doCaseManagement (pathogenesisModel->determineState (ageYears, withinHostModel),
				    withinHostModel,
				    latestReport,
				    ageYears,
				    _doomed);
}
