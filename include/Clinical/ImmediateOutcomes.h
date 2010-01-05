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

#ifndef Hmod_ClinicalImmediateOutcomes
#define Hmod_ClinicalImmediateOutcomes

#include "Clinical/ClinicalModel.h"
#include "OldCaseManagement.h"

namespace OM { namespace Clinical {

/** This implementation of the model is intended to use the old case-management
 * model with immediate outcomes of clinical events (immediate recovery with
 * total parasite clearance or immediate death). */
class ClinicalImmediateOutcomes : public ClinicalModel
{
public:
  /** Initialises parameters, loading from XML data. */
  static void initParameters ();
  
  ClinicalImmediateOutcomes (double cF, double tSF);
  ~ClinicalImmediateOutcomes ();
  
  void doClinicalUpdate (WithinHost::WithinHostModel& withinHostModel, double ageYears);
  
  inline bool recentTreatment() {
    return caseManagement->recentTreatment();
  }
  
protected:
  virtual void checkpoint (istream& stream);
  virtual void checkpoint (ostream& stream);
  
private:
  //TODO move implementation into class?
  /// The CaseManagementModel decides how to treat ill individuals
  OldCaseManagement * caseManagement;
};

} }
#endif
