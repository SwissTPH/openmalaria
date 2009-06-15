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


// -----  static init  -----

void ClinicalEventScheduler::init () {
  NewCaseManagement::init();
}


// -----  construction, destruction and checkpointing  -----

ClinicalEventScheduler::ClinicalEventScheduler (double cF, double tSF) :
    ClinicalModel (cF, tSF),
    caseManagement(new NewCaseManagement (tSF))
{}
ClinicalEventScheduler::~ClinicalEventScheduler() {
  delete caseManagement; 
}

ClinicalEventScheduler::ClinicalEventScheduler (istream& in) :
    ClinicalModel (in),
    caseManagement(new NewCaseManagement (in))
{}
void ClinicalEventScheduler::write (ostream& out) {
  pathogenesisModel->write (out);
  out << latestReport;
  out << _doomed << endl; 
  caseManagement->write (out);
}


// -----  other methods  -----

void ClinicalEventScheduler::doCaseManagement (WithinHostModel& withinHostModel, double ageYears) {
  caseManagement->doCaseManagement (pathogenesisModel->determineState (ageYears, withinHostModel),
				    withinHostModel,
				    latestReport,
				    ageYears,
				    _doomed);
}
