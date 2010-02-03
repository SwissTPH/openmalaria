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
#include "util/gsl.h"

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
    uint32_t start_bit = 0;
    
    const scnXml::DrugDescription& data = InputData.getScenario().getDrugDescription().get ();
    for (scnXml::DrugDescription::DrugConstIterator drug = data.getDrug().begin(); drug != data.getDrug().end(); ++drug) {
	DrugType::addDrug (new LSTMDrugType (*drug, start_bit));
    }
}

uint32_t LSTMDrugType::new_proteome_ID () {
    uint32_t id = 0;	// proteome / genotype identifier
    // for each drug / locus,
    for (map<string,DrugType*>::const_iterator it = available.begin(); it != available.end(); ++it) {
	const LSTMDrugType& dt = *(LSTMDrugType*)(it->second);
	double sample = gsl::rngUniform();
	for (size_t i = 0; i < dt.PD_params.size(); ++i) {
	    // we randomly pick an allele according to its initial frequency
	    if (sample <= dt.PD_params[i].cum_initial_frequency) {
		// add identifier for this allele into the proteome identifier
		assert ((~dt.allele_mask) & (i << dt.allele_rshift) == 0);	// sanity check (identifier within specified bits)
		id |= (i << dt.allele_rshift);
		goto gotAllele;	// and return
	    }
	}
	assert(false);	// sanity check that we did get an allele (only in debug mode)
	gotAllele:;
    }
    return id;	// done (includes components specifying each allele)
}



// -----  Non-static DrugType functions  -----

LSTMDrugType::LSTMDrugType (const scnXml::Drug& drugData, uint32_t& bit_start) :
    DrugType(drugData.getAbbrev()),
    allele_rshift (bit_start)
{
    const scnXml::PD::AlleleSequence& alleles = drugData.getPD().getAllele();
    if (alleles.size() < 1)
	throw util::xml_scenario_error ("Expected at least one allele for each drug.");
    
    // got length l; want minimal n such that: 2^n >= l
    // that is, n >= log_2 (l)
    // so n = ceil (log_2 (l))
    uint32_t n_bits = std::ceil (log (alleles.size()) / log(2.0));
    assert (std::pow (2, n_bits) >= alleles.size());
    allele_mask = uint32_t (std::pow (2, n_bits)) - 1;
    // update bit_start to next available bit:
    bit_start += n_bits;
    if (bit_start > 32)
	throw std::logic_error ("Implementation can't cope with this many alleles & drugs.");
    
    PD_params.resize (alleles.size());
    double cum_IF = 0.0;
    for (size_t i = 0; i < PD_params.size(); ++i) {
	cum_IF += alleles[i].getInitial_frequency ();
	PD_params[i].cum_initial_frequency = cum_IF;
	PD_params[i].max_killing_rate = alleles[i].getMax_killing_rate ();
	PD_params[i].IC50_pow_slope = pow(alleles[i].getIC50 (), alleles[i].getSlope ());
	PD_params[i].slope = alleles[i].getSlope ();
    }
    for (size_t i = 0; i < PD_params.size(); ++i) {
	PD_params[i].cum_initial_frequency /= cum_IF;	// scale: initial freq. of each is out of sum of all initial freq.s
    }
    // Be absolutely certain this is 1 so we can't get a random double greater than this value:
    PD_params[PD_params.size()-1].cum_initial_frequency = 1.0;
    
    elimination_rate_constant = log(2) / drugData.getPK().getHalf_life();
    vol_dist = drugData.getPK().getVol_dist();
}
LSTMDrugType::~LSTMDrugType () {}

} }