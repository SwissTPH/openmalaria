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

#include "PkPd/PkPdModel.h"
#include "PkPd/Drug/HoshenDrugType.h"
#include "Global.h"
#include "Host/proteome.h"

// submodels:
#include "PkPd/HoshenPkPdModel.h"
#include "PkPd/LSTMPkPdModel.h"


// Temporary switch to use the LSTM model − this may eventually be determined by XML data or XML model version.
const bool Use_LSTM = false;

// -----  static functions  -----

void PkPdModel::init () {
  if (Global::modelVersion & INCLUDES_PK_PD) {
    ProteomeManager::init ();
    if (Use_LSTM) {
      LSTMDrugType::init();
      LSTMDrug::init ();
      LSTMPkPdModel::init();
    }
    else {
      HoshenDrugType::init();
      HoshenDrug::init ();
      HoshenPkPdModel::init();
    }
  }
}
void PkPdModel::cleanup () {
    if (Global::modelVersion & INCLUDES_PK_PD) {
	ProteomeManager::cleanup ();
    }
}

void PkPdModel::readStatic (istream& in) {
  if (Global::modelVersion & INCLUDES_PK_PD) {
//     ProteomeManager::read (in);
  }
}
void PkPdModel::writeStatic (ostream& out) {
  if (Global::modelVersion & INCLUDES_PK_PD) {
//     ProteomeManager::write (out);
  }
}

PkPdModel* PkPdModel::createPkPdModel () {
  if (Global::modelVersion & INCLUDES_PK_PD) {
    if (Use_LSTM)
      return new LSTMPkPdModel ();
    else
      return new HoshenPkPdModel ();
  }
  return new PkPdModel();
}

PkPdModel* PkPdModel::createPkPdModel (istream& in) {
  if (Global::modelVersion & INCLUDES_PK_PD) {
    if (Use_LSTM)
      return new LSTMPkPdModel (in);
    else
      return new HoshenPkPdModel (in);
  }
  return new PkPdModel(in);
}

// -----  non-static functions  -----
