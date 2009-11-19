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

#ifndef Hmod_DUMMY_PK_PD_Drug_Model
#define Hmod_DUMMY_PK_PD_Drug_Model

#include "Drug/DrugModel.h"
#include "Drug/DummyPkPdDrug.h"

/** Pharmacokinetic and pharmacodynamics drug model. Dummy version
 *
 * (Excuse the horrible class name, but it allows keeping to the usual naming
 * convention.)
 * 
 * Holds per-human data for Tiago / the Liverpool school of medicine's
 * drug model.
 * 
 * Some of the implementation is contained in the drug.h/drug.cpp files. */
class DummyPkPdDrugModel : public DrugModel {
public:
  ///@brief Static functions
  //@{
  static void init ();
  //@}
  
  DummyPkPdDrugModel ();
  DummyPkPdDrugModel (istream& in);
  virtual ~DummyPkPdDrugModel ();
  virtual void write (ostream& out) const;
  
  void medicate(string drugAbbrev, double qty, int time, double age, double weight);
  void decayDrugs ();
  double getDrugFactor (const ProteomeInstance* infProteome);
  
private:
  list<DummyPkPdDrug> _drugs;
};

#endif
