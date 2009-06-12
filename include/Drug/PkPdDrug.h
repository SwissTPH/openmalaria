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

#ifndef Hmod_PK_PD_Drug
#define Hmod_PK_PD_Drug

#include "Drug/DrugModel.h"
#include "Drug/drug.h"

/** Pharmacokinetic and pharmacodynamics drug model.
 *
 * (Excuse the horrible class name, but it allows keeping to the usual naming
 * convention.)
 * 
 * Most of the implementation is contained in the drug.h/drug.cpp files, and
 * this class is just a wrapper to minimize changes to those files. */
class PkPdDrug : public DrugModel {
public:
  ///@brief Static functions
  //@{
  static void init ();
  //@}
  
  PkPdDrug ();
  PkPdDrug (istream& in);
  virtual ~PkPdDrug ();
  void write (ostream& out);
  
  void medicate(string drugAbbrev, double qty, int time);
  void decayDrugs ();
  
  void setWeight (double w);
  double getDrugFactor (ProteomeInstance* infProteome);
  
private:
  DrugProxy _proxy;
};

#endif
