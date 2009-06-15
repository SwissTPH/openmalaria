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

#ifndef Hmod_CaseManagementModel
#define Hmod_CaseManagementModel

#include "Pathogenesis/PathogenesisModel.h"
#include "event.h"

class WithinHostModel;

class CaseManagementModel
{
public:
  /// Report last event, if any
  virtual ~CaseManagementModel ();
  
  /** Determine treatment for a human.
   * @param pgState = Wellbeing of subject (well, severe malaria sickness, etc.)
   * @param withinHostModel = WithinHostModel of human.
   * @param ageYears = Age of human.
   * @param doomed = _doomed variable of Human; used to kill the human.
   *	Passing like this isn't ideal. */
  virtual void doCaseManagement (Pathogenesis::State pgState, WithinHostModel& withinHostModel, Event& latestReport, double ageYears, int& doomed) =0;
  
  bool recentTreatment();
  
  virtual void write(ostream& out) const;
  
  static int caseManagementMemory;
  
protected:
  /// Create, initializing _treatmentSeekingFactor to tSF.
  CaseManagementModel (double tSF);
  /// Create, loading from checkpoint
  CaseManagementModel (istream& in);
  
  //! treatment seeking for heterogeneity
  double _treatmentSeekingFactor;
  
  /** Timestep of the last treatment (TIMESTEP_NEVER if never treated). */
  int _tLastTreatment;
};

#endif
