/*

  This file is part of OpenMalaria.
 
  Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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

#include "PkPd/Drug/LSTMDrugType.h"
#include "inputData.h"
#include "util/errors.hpp"

#include <assert.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <sstream>

using namespace std;

namespace OM { namespace PkPd {
    
/*
 * Static variables and functions
 */

void LSTMDrugType::init () {
    DrugType::init();
    
    const scnXml::DrugDescription& data = InputData.getScenario().getDrugDescription().get ();
    for (scnXml::DrugDescription::DrugConstIterator drug = data.getDrug().begin(); drug != data.getDrug().end(); ++drug) {
	DrugType::addDrug (new LSTMDrugType (*drug));
    }
}

// -----  Non-static DrugType functions  -----

LSTMDrugType::LSTMDrugType (const scnXml::Drug& drugData)
: DrugType(drugData.getAbbrev())
{
    const scnXml::PD::AlleleSequence& alleles = drugData.getPD().getAllele();
    if (alleles.size() < 1)
	throw util::xml_scenario_error ("Expected at least one allele for each drug.");
    PD_params.resize (alleles.size());
    for (size_t i = 0; i < alleles.size(); ++i) {
	PD_params[i].max_killing_rate = alleles[i].getMax_killing_rate ();
	PD_params[i].IC50 = alleles[i].getIC50 ();
	PD_params[i].slope = alleles[i].getSlope ();
    }
    
    elimination_rate_constant = log(2) / drugData.getPK().getHalf_life();
    vol_dist = drugData.getPK().getVol_dist();
}
LSTMDrugType::~LSTMDrugType () {}

} }