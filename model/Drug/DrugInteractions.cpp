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

#include "Drug/DrugInteractions.h"
#include "Global.h"
#include "proteome.h"

// submodels:
#include "Drug/DummyPkPdDrugInteractions.h"


// -----  static functions  -----

void DrugInteractions::init () {
  if (Global::modelVersion & INCLUDES_PK_PD) {
    initProteomeModule();
    DummyPkPdDrugInteractions::init();
  }
}

void DrugInteractions::readStatic (istream& in) {
  if (Global::modelVersion & INCLUDES_PK_PD) {
    ProteomeManager::read (in);
  }
}
void DrugInteractions::writeStatic (ostream& out) {
  if (Global::modelVersion & INCLUDES_PK_PD) {
    ProteomeManager::write (out);
  }
}

DrugInteractions* DrugInteractions::createDrugInteractions () {
  if (Global::modelVersion & INCLUDES_PK_PD) {
    return new DummyPkPdDrugInteractions ();
  }
  return new DrugInteractions();
}

DrugInteractions* DrugInteractions::createDrugInteractions (istream& in) {
  if (Global::modelVersion & INCLUDES_PK_PD) {
    return new DummyPkPdDrugInteractions (in);
  }
  return new DrugInteractions(in);
}

// -----  non-static functions  -----
