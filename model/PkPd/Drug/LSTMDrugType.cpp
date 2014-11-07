/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2014 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2014 Liverpool School Of Tropical Medicine
 * 
 * OpenMalaria is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "PkPd/Drug/LSTMDrugType.h"
#include "WithinHost/Genotypes.h"
#include "util/errors.h"
#include "util/random.h"
#include <schema/pharmacology.h>

#include <cmath>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/format.hpp>

using namespace std;

namespace OM {
namespace PkPd {
using namespace OM::util;
using WithinHost::Genotypes;

// -----  Static variables and functions  -----

// The list of drugTypes drugs. Not checkpointed.
//TODO: do we need to store by pointer?
typedef boost::ptr_vector<LSTMDrugType> DrugTypesT;
DrugTypesT drugTypes;
map<string,size_t> drugTypeNames;


void LSTMDrugType::init (const scnXml::Pharmacology::DrugsType& drugData) {
    foreach( const scnXml::PKPDDrug& drug, drugData.getDrug() ){
        const string& abbrev = drug.getAbbrev();
        // Check drug doesn't already exist
        if (drugTypeNames.find (abbrev) != drugTypeNames.end())
            throw TRACED_EXCEPTION_DEFAULT (string ("Drug added twice: ").append(abbrev));
        
        size_t i = drugTypes.size();
        drugTypes.push_back( new LSTMDrugType (i, drug) );
        drugTypeNames[abbrev] = i;
    }
}
void LSTMDrugType::clear()
{
    drugTypes.clear();
    drugTypeNames.clear();
}

size_t LSTMDrugType::findDrug(string _abbreviation) {
    map<string,size_t>::const_iterator it = drugTypeNames.find (_abbreviation);
    if (it == drugTypeNames.end())
        throw util::xml_scenario_error (string ("attempt to use drug without description: ").append(_abbreviation));
    
    return it->second;
}

const LSTMDrugType& LSTMDrugType::getDrug(size_t index) {
    return drugTypes[index];
}



// -----  Non-static LSTMDrugType functions  -----

LSTMDrugType::LSTMDrugType (size_t index, const scnXml::PKPDDrug& drugData) :
        index (index),
        name(drugData.getAbbrev())
{
    const scnXml::PD::PhenotypeSequence& pElt = drugData.getPD().getPhenotype();
    assert( pElt.size() > 0 );  // required by XSD
    
    negligible_concentration = drugData.getPK().getNegligible_concentration();
    neg_elimination_rate_constant = -log(2.0) / drugData.getPK().getHalf_life();
    vol_dist = drugData.getPK().getVol_dist();
    
    set<string> loci_per_phenotype;
    string first_phenotype_name;
    size_t n_phenotypes = pElt.size();
    PD.reserve (n_phenotypes);
    // per phenotype (first index), per locus-of-restriction (second list), a list of alleles
    vector<vector<vector<uint32_t> > > phenotype_restrictions;
    phenotype_restrictions.reserve( n_phenotypes );
    for (size_t i = 0; i < n_phenotypes; ++i) {
        const scnXml::Phenotype& phenotype = pElt[i];
        if( i == 0 ){
            if( phenotype.getName().present() )
                first_phenotype_name = phenotype.getName().get();
            else
                first_phenotype_name = (boost::format("(%1%)") %i).str();
            foreach( const scnXml::Restriction& restriction, phenotype.getRestriction() ){
                loci_per_phenotype.insert( restriction.getOnLocus() );
            }
        }else{
            set<string> expected_loci = loci_per_phenotype;
            foreach( const scnXml::Restriction& restriction, phenotype.getRestriction() ){
                size_t n = expected_loci.erase( restriction.getOnLocus() );
                if( n == 0 ){
                    string this_phenotype_name;
                    if( phenotype.getName().present() )
                        this_phenotype_name = phenotype.getName().get();
                    else
                        this_phenotype_name = (boost::format("(%1%)") %i).str();
                    throw util::xml_scenario_error( (boost::format(
                        "pharmacology/drugs/drug/pd/phenotype/restriction/onLocus:"
                        "locus %1% included in restriction of phenotype \"%2%\" but not %3%")
                        %restriction.getOnLocus() %this_phenotype_name
                        %first_phenotype_name).str() );
                }
            }
            if( !expected_loci.empty() ){
                string this_phenotype_name;
                if( phenotype.getName().present() )
                    this_phenotype_name = phenotype.getName().get();
                else
                    this_phenotype_name = (boost::format("(%1%)") %i).str();
                throw util::xml_scenario_error( (boost::format(
                    "pharmacology/drugs/drug/pd/phenotype/restriction/onLocus:"
                    "locus %1% included in restriction of phenotype \"%2%\" but not %3%")
                    %(*expected_loci.begin()) %first_phenotype_name
                    %this_phenotype_name).str() );
            }
        }
        map<string,size_t> loci;
        vector<vector<uint32_t> > loc_alleles;
        foreach( const scnXml::Restriction& restriction, phenotype.getRestriction() ){
            uint32_t allele = Genotypes::findAlleleCode( restriction.getOnLocus(), restriction.getToAllele() );
            if( allele == numeric_limits<uint32_t>::max() ){
                throw util::xml_scenario_error( (boost::format("phenotype has "
                    "restriction on locus %1%, allele %2% but this locus/allele "
                    "has not been defined in parasiteGenetics section")
                    %restriction.getOnLocus() %restriction.getToAllele()).str() );
            }
            map<string,size_t>::const_iterator it = loci.find(restriction.getOnLocus());
            if( it == loci.end() ){
                loci[restriction.getOnLocus()] = loc_alleles.size();
                loc_alleles.push_back( vector<uint32_t>(1,allele) );
            }else{
                loc_alleles[it->second].push_back( allele );
            }
        }
        phenotype_restrictions.push_back( loc_alleles );
        PD.push_back( new LSTMDrugPD( pElt[i], -neg_elimination_rate_constant ) );
    }
    
    if( loci_per_phenotype.size() == 0 ){
        if( pElt.size() > 1 ){
            throw util::xml_scenario_error( "pharmacology/drugs/drug/pd/phenotype:"
                " restrictions required when num. phenotypes > 1" );
        }else{
            // All genotypes map to phenotype 0
            genotype_mapping.assign( Genotypes::getGenotypes().size(), 0 );
        }
    }else{
        const vector<Genotypes::Genotype>& genotypes = Genotypes::getGenotypes();
        genotype_mapping.assign( genotypes.size(), 0 );
        for( size_t j = 0; j < genotypes.size(); ++j ){
            uint32_t phenotype = numeric_limits<uint32_t>::max();
            for( size_t i = 0; i < n_phenotypes; ++i ){
                // genotype matches this phenotype when, for every locus among
                // the phenotype restriction rules, one allele matches a
                // genotype allele
                bool match = true;
                foreach( const vector<uint32_t>& loc_alleles, phenotype_restrictions[i] ){
                    bool match_for_locus = false;
                    foreach( uint32_t restrict_allele, loc_alleles ){
                        if( genotypes[j].alleles.count( restrict_allele ) > 0 ){
                            assert( !match_for_locus ); // shouldn't be two matches
                            match_for_locus = true;
                        }
                    }
                    if( !match_for_locus ){
                        match = false;
                        break;  // can skip rest of rules
                    }
                }
                if( match ){
                    if( phenotype == numeric_limits<uint32_t>::max() ){
                        phenotype = i;
                    }else{
                        // We could try to convert genotype's allele codes to
                        // locus/allele names, but requires several lookups
                        throw util::xml_scenario_error("phenotype restrictions "
                            "not restrictive enough: multiple phenotypes match "
                            "some genotype");
                    }
                }
            }
            if( phenotype == numeric_limits<uint32_t>::max() ){
                // We could try to convert genotype's allele codes to
                // locus/allele names, but requires several lookups
                throw util::xml_scenario_error("phenotype restrictions too"
                    "restrictive: no phenotype matching some genotype");
            }
            genotype_mapping[j] = phenotype;
        }
    }
}
LSTMDrugType::~LSTMDrugType () {
}

const LSTMDrugPD& LSTMDrugType::getPD( uint32_t genotype ) const {
    return PD[genotype_mapping[genotype]];
}

void LSTMDrugType::updateConcentration( double& C0, double duration ) const {
    // exponential decay of drug concentration
    C0 *= exp(neg_elimination_rate_constant * duration);
}
void LSTMDrugType::updateConcentrationIV( double& C0, double duration, double rate ) const {
    // exponential decay of drug concentration
    C0 *= exp(neg_elimination_rate_constant * duration);
    // TODO: explain this
    // TODO: why not adjust rate by neg_elimination_rate_constant and/or vol_dist earlier?
    C0 += rate
          * (1.0 - exp(neg_elimination_rate_constant * duration) )
          / ( -neg_elimination_rate_constant * vol_dist );
}

}
}
