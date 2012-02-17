/*
  This file is part of OpenMalaria.

  Copyright (C) 2005-2011 Swiss Tropical Institute and Liverpool School Of Tropical Medicine

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
#include "util/errors.h"
#include "util/random.h"

#include <cmath>

using namespace std;

namespace OM {
namespace PkPd {
using namespace OM::util;

// -----  Static variables and functions  -----

LSTMDrugType::Available LSTMDrugType::available;


void LSTMDrugType::init (const scnXml::Pharmacology& data) {
    uint32_t start_bit = 0;

    for (scnXml::Pharmacology::DrugConstIterator drug = data.getDrug().begin(); drug != data.getDrug().end(); ++drug) {
        LSTMDrugType::addDrug (auto_ptr<LSTMDrugType>(new LSTMDrugType (*drug, start_bit)));
    }
}
void LSTMDrugType::cleanup () {
    for ( Available::const_iterator it = available.begin(); it != available.end(); ++it ) {
        delete it->second;
    }
    available.clear ();
}

void LSTMDrugType::addDrug(auto_ptr<LSTMDrugType> drug) {
    const string& abbrev = drug->abbreviation;
    // Check drug doesn't already exist
    if (available.find (abbrev) != available.end())
        throw TRACED_EXCEPTION_DEFAULT (string ("Drug already in registry: ").append(abbrev));

    available[ abbrev ] = drug.release();
}

const LSTMDrugType& LSTMDrugType::getDrug(string _abbreviation) {
    Available::const_iterator i = available.find (_abbreviation);
    if (i == available.end())
        throw util::xml_scenario_error (string ("prescribed non-existant drug ").append(_abbreviation));

    return *i->second;
}

uint32_t LSTMDrugType::new_proteome_ID () {
    uint32_t id = 0;    // proteome / genotype identifier
    // for each drug / locus,
    for (Available::const_iterator it = available.begin(); it != available.end(); ++it) {
        const LSTMDrugType& dt = *it->second;
        double sample = random::uniform_01();
        for (size_t i = 0; i < dt.drugAllele.size(); ++i) {
            // we randomly pick an allele according to its initial frequency
            if (sample <= dt.cumInitialFreq[i]) {
                // add identifier for this allele into the proteome identifier
                assert (((~dt.allele_mask) & (i << dt.allele_rshift)) == 0);    // sanity check (identifier within specified bits)
                id |= (i << dt.allele_rshift);
                goto gotAllele; // and return
            }
        }
        assert(false);  // sanity check that we did get an allele (only in debug mode)
        gotAllele: ;
    }
    return id;  // done (includes components specifying each allele)
}



// -----  Non-static LSTMDrugType functions  -----

LSTMDrugType::LSTMDrugType (const scnXml::Drug& drugData, uint32_t& bit_start) :
        abbreviation (drugData.getAbbrev()),
        allele_rshift (bit_start)
{
    const scnXml::PD::AlleleSequence& alleles = drugData.getPD().getAllele();
    if (alleles.size() < 1)
        throw util::xml_scenario_error ("Expected at least one allele for each drug.");

    // got length l; want minimal n such that: 2^n >= l
    // that is, n >= log_2 (l)
    // so n = ceil (log_2 (l))
    uint32_t n_bits = (uint32_t)std::ceil (log (double(alleles.size())) / log(2.0));
    assert (std::pow (2.0, (double)n_bits) >= alleles.size());
    allele_mask = static_cast<uint32_t>((std::pow (2.0, (double)n_bits)) - 1);
    // update bit_start to next available bit:
    bit_start += n_bits;
    if (bit_start > 32)
        throw TRACED_EXCEPTION_DEFAULT ("Implementation can't cope with this many alleles & drugs.");

    negligible_concentration = drugData.getPK().getNegligible_concentration();
    neg_elimination_rate_constant = -log(2.0) / drugData.getPK().getHalf_life();
    vol_dist = drugData.getPK().getVol_dist();

    cumInitialFreq.reserve(alleles.size());
    double cum_IF = 0.0;
    for (size_t i = 0; i < alleles.size(); ++i) {
        cum_IF += alleles[i].getInitial_frequency ();
        cumInitialFreq.push_back( cum_IF );
    }
    for (size_t i = 0; i < cumInitialFreq.size(); ++i) {
        // scale: initial freq. of each is out of sum of all initial freq.s
        cumInitialFreq[i] /= cum_IF;
    }
    // Be absolutely certain this is 1 so we can't get a random double greater than this value:
    cumInitialFreq[cumInitialFreq.size()-1] = 1.0;

    drugAllele.reserve (alleles.size());
    for (size_t i = 0; i < alleles.size(); ++i) {
        drugAllele.push_back( new LSTMDrugAllele( alleles[i], -neg_elimination_rate_constant ) );
    }
}
LSTMDrugType::~LSTMDrugType () {
}

const LSTMDrugAllele& LSTMDrugType::getAllele( uint32_t proteome_ID ) const {
    uint32_t allele = (proteome_ID >> allele_rshift) & allele_mask;
    return drugAllele[allele];
}

void LSTMDrugType::updateConcentration( double& C0, double duration ) const {
    C0 *= exp(neg_elimination_rate_constant * duration);
}
void LSTMDrugType::updateConcentrationIV( double& C0, double duration, double rate ) const {
    C0 *= exp(neg_elimination_rate_constant * duration);
    C0 += rate
          * (1.0 - exp(neg_elimination_rate_constant * duration) )
          / ( -neg_elimination_rate_constant * vol_dist );
}

}
}
