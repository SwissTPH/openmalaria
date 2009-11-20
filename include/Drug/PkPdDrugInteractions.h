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

#ifndef Hmod_PK_PD_Drug_Interactions
#define Hmod_PK_PD_Drug_Interactions

#include "Drug/DrugInteractions.h"
#include "Drug/PkPdDrug.h"

/** Pharmacokinetic and pharmacodynamics drug Interactions. 
 *
 * 
 * Some of the implementation is contained in the drug.h/drug.cpp files. */
class PkPdDrugInteractions : public DrugInteractions {
public:
  ///@brief Static functions
  //@{
  static void init ();
  //@}
  
  PkPdDrugInteractions ();
  PkPdDrugInteractions (istream& in);
  virtual ~PkPdDrugInteractions ();
  virtual void write (ostream& out) const;
  
  void medicate(string drugAbbrev, double qty, int time, double age, double weight);
  void decayDrugs ();
  double getDrugFactor (const ProteomeInstance* infProteome);
  
private:
  list<PkPdDrug> _drugs;
};

#endif
