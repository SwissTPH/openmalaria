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

#include "morbidityModel.h"
#include "event.h"

class WithinHostModel;

class CaseManagementModel
{
public:
  /// Initialize caseManagementMemory
  static void init();
  /// Return a new CaseManagementModel. Don't create one directly.
  static CaseManagementModel* createCaseManagementModel (double tSF);

  /// Report last event, if any
  ~CaseManagementModel ();
  
  /** Determine treatment for a human.
   * @param infection = Type of infection
   * @param withinHostModel = WithinHostModel of human.
   * @param ageYears = Age of human.
   * @param doomed = _doomed variable of human. Passing like this isn't ideal.
   */
  virtual void doCaseManagement (Morbidity::Infection infection, WithinHostModel& withinHostModel, double ageYears, int& doomed) =0;
  
  bool recentTreatment();
  
  /** Return the case management's event.
   * NOTE: possibly this method should be removed later. */
  Event& getEvent();
  
  virtual void write(ostream& out) const;
  virtual void read(istream& in);
  
  static int caseManagementMemory;
  
protected:
  /// Create, initializing _treatmentSeekingFactor to tSF.
  CaseManagementModel (double tSF);
  
  /** Next event to report.
   * Only reported when the Human dies or a separate episode occurs. */
  Event _latestEvent;
  //! treatment seeking for heterogeneity
  double _treatmentSeekingFactor;
  
  //!time of the last treatment
  int _tLastTreatment;
};

#endif
