/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2015 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
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
#include "util/CommandLine.h"
#include "PkPd/Drug/LSTMDrugOneComp.h"
#include "PkPd/Drug/LSTMDrugThreeComp.h"
#include "PkPd/Drug/LSTMDrugConversion.h"
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
// List of all indices of drugs being used
vector<size_t> drugsInUse;


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

size_t LSTMDrugType::numDrugTypes(){
    return drugTypes.size();
}

// Add index to drugsInUse if not already present. Not fast but doesn't need to be.
void drugIsUsed(size_t index){
    foreach( size_t i, drugsInUse ){
        if( i == index ) return;        // already in list
    }
    drugsInUse.push_back(index);
}
size_t LSTMDrugType::findDrug(string _abbreviation) {
    map<string,size_t>::const_iterator it = drugTypeNames.find (_abbreviation);
    if (it == drugTypeNames.end())
        throw util::xml_scenario_error (string ("attempt to use drug without description: ").append(_abbreviation));
    size_t index = it->second;
    
    // We assume that drugs are used when and only when findDrug returns their
    // index or they are a metabolite of a drug returned here.
    drugIsUsed(index);
    if( drugTypes[index].conversion_rate.isSet() /*conversion model used*/ ){
        drugIsUsed(drugTypes[index].metabolite);
    }
    
    return index;
}

const vector< size_t >& LSTMDrugType::getDrugsInUse(){
    return drugsInUse;
}

LSTMDrug* LSTMDrugType::createInstance(size_t index) {
    LSTMDrugType& typeData = drugTypes[index];
    if( typeData.conversion_rate.isSet() ){
        LSTMDrugType& metaboliteData = drugTypes[typeData.metabolite];
        return new LSTMDrugConversion( typeData, metaboliteData );
    }else if( typeData.a12.isSet() ){
        // a21 is set when a12 is set; a13 and a31 may be set
        return new LSTMDrugThreeComp( typeData );
    }else{
        // none of a12/a21/a13/a31 should be set in this case
        return new LSTMDrugOneComp( typeData );
    }
}



// -----  Non-static LSTMDrugType functions  -----

LSTMDrugType::LSTMDrugType (size_t index, const scnXml::PKPDDrug& drugData) :
        index (index), metabolite(0),
        negligible_concentration(numeric_limits<double>::quiet_NaN()),
        neg_m_exp(numeric_limits<double>::quiet_NaN()),
        mwr(numeric_limits<double>::quiet_NaN())
{
    // ———  PK parameters  ———
    const scnXml::PK& pk = drugData.getPK();
    negligible_concentration = pk.getNegligible_concentration();
    if( pk.getHalf_life().present() ){
        elimination_rate.setParams( log(2.0) / pk.getHalf_life().get(), 0.0 );
        neg_m_exp = 0.0; // no dependence on body mass
    }else{
        if( !(pk.getK().present() && pk.getM_exponent().present()) ){
            throw util::xml_scenario_error( "PK data must include either half_life or (k and m_exponent)" );
        }
        elimination_rate.setParams( pk.getK().get() );
        neg_m_exp = -pk.getM_exponent().get();
    }
    vol_dist.setParams( pk.getVol_dist(),
                        pk.getVol_dist().getSigma() );
    if( pk.getCompartment2().present() ){
        a12.setParams( pk.getCompartment2().get().getA12() );
        a21.setParams( pk.getCompartment2().get().getA21() );
        if( pk.getCompartment3().present() ){
            a13.setParams( pk.getCompartment3().get().getA13() );
            a31.setParams( pk.getCompartment3().get().getA31() );
        }else{
            // 2-compartment model: use 3-compartment code with these parameters set to zero
            a13.setParams(0.0, 0.0);
            a31.setParams(0.0, 0.0);
        }
    }else if( pk.getCompartment3().present() ){
        throw util::xml_scenario_error( "PK model specifies parameters for "
                "compartment3 without compartment2" );
    }
    if( pk.getConversion().present() ){
        if( a12.isSet() ){
            throw util::xml_scenario_error( "PK conversion model is incompatible with 2/3-compartment model" );
        }
        const scnXml::Conversion& conv = pk.getConversion().get();
        try{
            metabolite = findDrug(conv.getMetabolite());
        }catch( util::xml_scenario_error e ){
            throw util::xml_scenario_error( "PK: metabolite drug not found; metabolite must be defined *before* parent drug!" );
        }
        conversion_rate.setParams( conv.getRate() );
        mwr = conv.getMolRatio();
    }
    if( pk.getK_a().present() ){
        if( !a12.isSet() && !conversion_rate.isSet() ){
            throw util::xml_scenario_error( "PK models only allow an "
                "absorption rate parameter (k_a) when compartment2 or "
                " conversion parameters are present" );
        }
        absorption_rate.setParams(pk.getK_a().get());
    }else{
        if( a12.isSet() ){
            throw util::xml_scenario_error( "PK models require an absorption "
                "rate parameter (k_a) when compartment2 is present" );
        }
    }
    
    // ———  PD parameters  ———
    const scnXml::PD::PhenotypeSequence& pd = drugData.getPD().getPhenotype();
    assert( pd.size() > 0 );  // required by XSD
    
    set<string> loci_per_phenotype;
    string first_phenotype_name;
    size_t n_phenotypes = pd.size();
    PD.reserve (n_phenotypes);
    // per phenotype (first index), per locus-of-restriction (second list), a list of alleles
    vector<vector<vector<uint32_t> > > phenotype_restrictions;
    phenotype_restrictions.reserve( n_phenotypes );
    for(size_t i = 0; i < n_phenotypes; ++i) {
        const scnXml::Phenotype& phenotype = pd[i];
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
        PD.push_back( new LSTMDrugPD( pd[i] ) );
    }
    
    if( loci_per_phenotype.size() == 0 ){
        if( pd.size() > 1 ){
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
        #ifdef WITHOUT_BOINC
        if( util::CommandLine::option( util::CommandLine::PRINT_GENOTYPES ) ){
            cout << endl;
            // Debug output to see which genotypes correspond to phenotypes
            cout << "Phenotype mapping: " << endl << "----------" << endl;
            boost::format fmtr("|%8d|%14d|");
            cout << (fmtr % "genotype" % "phenotype") << endl;
            cout << (fmtr % "--------" % "-------------") << endl;
            uint32_t genotype = 0;
            stringstream phenotypeName;
            for( vector<uint32_t>::const_iterator phen = genotype_mapping.begin(); phen !=genotype_mapping.end(); ++phen ){
                phenotypeName.str("");
                if(pd[*phen].getName().present()){
                    phenotypeName << pd[*phen].getName().get();
                } else {
                    phenotypeName << "no: " << *phen;
                }
                cout << (fmtr % (genotype*100000) % phenotypeName.str() ) << endl;
                genotype += 1;
            }
        #endif
        }
    }
}
LSTMDrugType::~LSTMDrugType () {
}

const LSTMDrugPD& LSTMDrugType::getPD( uint32_t genotype ) const {
    return PD[genotype_mapping[genotype]];
}

}
}
