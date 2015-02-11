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

#include "mon/info.h"
#include "mon/reporting.h"
#include "mon/management.h"
#define H_OM_mon_cpp
#include "mon/OutputMeasures.hpp"
#include "WithinHost/Genotypes.h"
#include "Clinical/CaseManagementCommon.h"
#include "Host/Human.h"
#include "util/errors.h"
#include "schema/scenario.h"

#include <FastDelegate.h>
#include <iostream>
#include <boost/format.hpp>

namespace OM {
namespace mon {
    using namespace fastdelegate;

namespace impl {
    // Accumulators, variables:
    size_t currentSurvey = NOT_USED;
    SimTime nextSurveyTime = sim::future();
}

#ifndef NDEBUG
const size_t NOT_ACCEPTED = boost::integer_traits<size_t>::const_max - 1;
#endif

typedef FastDelegate4<ostream&,size_t,int,size_t> WriteDelegate;

// Store by human age and cohort
template<typename T, bool BY_AGE, bool BY_COHORT, bool BY_SPECIES,
    bool BY_GENOTYPE, bool BY_DRUG>
class Store{
    // This maps from an index in reports to an output measure
    vector<int> outMeasures;
    
    // mIndices and deployIndices both map from measures to indices.
    // The former should be faster but is insufficient for deployments.
    // Usage should not overlap (i.e. only one should match any measure).
    
    // This maps from measures (MHR_HOSTS, etc.) to an index in reports or NOT_USED
    vector<size_t> mIndices;
    // This maps from measures (MHD_VACCINATIONS, etc.) to a Deploy::Method enum
    // and an index in reports. Measures may have any number of matches here.
    typedef multimap<int,pair<uint8_t,size_t> > DeployInd_t;
    DeployInd_t deployIndices;
    
    size_t nAgeGroups, nCohortSets, nSpecies, nGenotypes, nDrugs;
    // These are the stored reports (multidimensional; use size() and index())
    vector<T> reports;
    
    // get size of reports
    inline size_t size(){ return outMeasures.size() * impl::nSurveys *
        nAgeGroups * nCohortSets * nSpecies * nGenotypes * nDrugs; }
    // get an index in reports
    inline size_t index( size_t m, size_t s, size_t a, size_t c, size_t sp, size_t g, size_t d ){
        return d + nDrugs *
            (g + nGenotypes *
            (sp + nSpecies *
            (c + nCohortSets *
            (a + nAgeGroups *
            (s + impl::nSurveys *
            m)))));
    }
    
    inline void add(T val, size_t mIndex, size_t survey, size_t ageIndex,
                    uint32_t cohortSet, size_t species, size_t genotype, size_t drug)
    {
#ifndef NDEBUG
        if( mIndex >= outMeasures.size() ||
            survey >= impl::nSurveys ||
            ageIndex >= nAgeGroups ||
            cohortSet >= nCohortSets ||
            species >= nSpecies ||
            genotype >= nGenotypes ||
            drug >= nDrugs
        ){
            cout << "Index out of bounds for survey:\t" << survey << " of " << impl::nSurveys
                << "\nmeasure number\t" << mIndex << " of " << outMeasures.size()
                << "\nage group\t" << ageIndex << " of " << nAgeGroups
                << "\ncohort set\t" << cohortSet << " of " << nCohortSets
                << "\nspecies\t" << species << " of " << nSpecies
                << "\ngenotype\t" << genotype << " of " << nGenotypes
                << "\ndrug\t" << drug << " of " << nDrugs
                << endl;
        }
#endif
        reports[index(mIndex,survey,ageIndex,cohortSet,species,genotype,drug)] += val;
    }
    
    void writeM( ostream& stream, size_t survey, int outMeasure, size_t inMeasure ){
        assert( !(BY_DRUG && (BY_SPECIES || BY_GENOTYPE)) );
        if( BY_SPECIES ){
            assert( !BY_AGE && !BY_COHORT );  // output col2 conflicts
            for( size_t species = 0; species < nSpecies; ++species ){
            for( size_t genotype = 0; genotype < nGenotypes; ++genotype ){
                const int col2 = species + 1 +
                    1000000 * genotype;
                T value = reports[index(inMeasure,survey,0,0,species,genotype,0)];
                stream << (survey+1) << '\t' << col2 << '\t' << outMeasure
                    << '\t' << value << lineEnd;
            } }
        }else if( BY_DRUG ){
            // Backwards compatibility: first age group starts at 1, unless
            // there isn't an age group:
            const int ageGroupAdd = BY_AGE ? 1 : 0;
            for( size_t cohortSet = 0; cohortSet < nCohortSets; ++cohortSet ){
            for( size_t ageGroup = 0; ageGroup < nAgeGroups; ++ageGroup ){
            for( size_t drug = 0; drug < nDrugs; ++drug ){
                // Yeah, >999 age groups clashes with cohort sets, but unlikely a real issue
                const int col2 = ageGroup + ageGroupAdd +
                    1000 * internal::cohortSetOutputId( cohortSet ) +
                    1000000 * (drug + 1);
                T value = reports[index(inMeasure,survey,ageGroup,cohortSet,0,0,drug)];
                stream << (survey+1) << '\t' << col2 << '\t' << outMeasure
                    << '\t' << value << lineEnd;
            } } }
        }else{
            // Backwards compatibility: first age group starts at 1, unless
            // there isn't an age group:
            const int ageGroupAdd = BY_AGE ? 1 : 0;
            for( size_t cohortSet = 0; cohortSet < nCohortSets; ++cohortSet ){
            for( size_t ageGroup = 0; ageGroup < nAgeGroups; ++ageGroup ){
            for( size_t genotype = 0; genotype < nGenotypes; ++genotype ){
                // Yeah, >999 age groups clashes with cohort sets, but unlikely a real issue
                const int col2 = ageGroup + ageGroupAdd +
                    1000 * internal::cohortSetOutputId( cohortSet ) +
                    1000000 * genotype;
                T value = reports[index(inMeasure,survey,ageGroup,cohortSet,0,genotype,0)];
                stream << (survey+1) << '\t' << col2 << '\t' << outMeasure
                    << '\t' << value << lineEnd;
            } } }
        }
    }
    
public:
    // Set up ready to accept reports.
    void init( const list<OutMeasure>& required, size_t nSp, size_t nD ){
        // Set these internally to 1 if we don't segregate.
        // Age groups -1 because last isn't reported.
        nAgeGroups = BY_AGE ? (AgeGroup::numGroups() - 1) : 1;
        nCohortSets = BY_COHORT ? impl::nCohortSets : 1;
        nSpecies = BY_SPECIES ? nSp : 1;
        nGenotypes = BY_GENOTYPE ? WithinHost::Genotypes::N() : 1;
        nDrugs = BY_DRUG ? nD : 1;
        mIndices.assign( M_NUM, NOT_USED );
        // outMeasures.size() is the number of measures we store here
        outMeasures.assign( 0, 0 );
        
        for( list<OutMeasure>::const_iterator it = required.begin();
            it != required.end(); ++it )
        {
            if( it->m >= M_NUM ) continue;      // skip: obsolete/special
            if( it->isDouble != (typeid(T) == typeid(double) ) ){
#ifndef NDEBUG
                // Debug mode: this should prevent silly errors where the type
                // reported does not match the type defined for some output:
                mIndices[it->m] = NOT_ACCEPTED;
#endif
                continue;
            }
            if( it->byAge != BY_AGE ||
                it->byCohort != BY_COHORT ||
                it->bySpecies != BY_SPECIES ||
                it->byGenotype != BY_GENOTYPE ||
                it->byDrug != BY_DRUG ) continue;
            
            if( mIndices[it->m] != NOT_USED ||
                (it->method == Deploy::NA && deployIndices.count(it->m) > 0) )
            {
                //NOTE: if we give MHR_HOSTS, etc. names visible to the XML
                // we should report that name. Current use of a number may be confusing.
                ostringstream msg;
                msg << "multiple use of monitoring measure " << it->m << " (used by ";
                findNamedMeasuresUsing( it->m, msg );
                msg << ") by age and cohort";
                throw util::xml_scenario_error( msg.str() );
            }
            
            size_t newInd = outMeasures.size();     // length becomes our index
            if( it->method == Deploy::NA ){
                mIndices[it->m] = newInd;
            }else{
                deployIndices.insert( make_pair(it->m, make_pair(it->method, newInd)) );
            }
            outMeasures.push_back( it->outId ); // increment length
        }
        
        reports.assign(size(), 0);
    }
    
    // Take a reported value and either store it or forget it.
    // If some of ageIndex, cohortSet, species are not applicable, use 0.
    void report( T val, Measure measure, size_t survey, size_t ageIndex,
                 uint32_t cohortSet, size_t species, size_t genotype, size_t drug )
    {
        if( survey == NOT_USED ) return; // pre-main-sim & unit tests we ignore all reports
        // last category is for humans too old for reporting groups:
        if( ageIndex == nAgeGroups ) return;
        assert( measure < mIndices.size() );
        assert( mIndices[measure] != NOT_ACCEPTED );
        assert( ageIndex < nAgeGroups && (BY_AGE || nAgeGroups == 1) );
        assert( cohortSet < nCohortSets && (BY_COHORT || nCohortSets == 1) );
        assert( species < nSpecies && (BY_SPECIES || nSpecies == 1) );
        assert( genotype < nGenotypes && (BY_GENOTYPE || nGenotypes == 1) );
        assert( drug < nDrugs && (BY_DRUG || nDrugs == 1) );
        if( mIndices[measure] == NOT_USED ){    // measure not used by this store
            assert( deployIndices.count(measure) == 0 );
            return;
        }
        add( val, mIndices[measure], survey, ageIndex, cohortSet, species,
             genotype, drug );
    }
    
    // Take a deployment report and potentially store it in one or more places
    // If some ageIndex or cohortSet are not applicable, use 0.
    void deploy( T val, Measure measure, size_t survey, size_t ageIndex,
                 uint32_t cohortSet, Deploy::Method method )
    {
        if( survey == NOT_USED ) return; // pre-main-sim & unit tests we ignore all reports
        if( ageIndex == nAgeGroups ) return;    // last category is for humans too old for reporting groups
        assert( method == Deploy::NA || method == Deploy::TIMED ||
            method == Deploy::CTS || method == Deploy::TREAT );
        assert( measure < mIndices.size() );
        assert( mIndices[measure] == NOT_USED );
        assert( nSpecies == 1 && nGenotypes == 1 );     // never used for deployments
        
        typedef DeployInd_t::const_iterator const_it_t;
        pair<const_it_t,const_it_t> range = deployIndices.equal_range( measure );
        for( const_it_t it = range.first; it != range.second; ++it ){
            uint8_t mask = it->second.first;
            if( (mask & method) == 0 ) continue;
            add( val, it->second.second, survey, ageIndex, cohortSet, 0, 0, 0 );
        }
    }
    
    // Return true if reports by this measure are recorded, false if they are discarded.
    bool isUsed( Measure measure ){
        assert( measure < mIndices.size() );
        assert( mIndices[measure] != NOT_ACCEPTED );
        return mIndices[measure] != NOT_USED || deployIndices.count(measure) > 0;
    }
    
    // Order self in a list of outputs
    void addMeasures( map<int,pair<WriteDelegate,size_t> >& mOrdered ){
        for( size_t m = 0; m < outMeasures.size(); ++m ){
            mOrdered[outMeasures[m]] = make_pair( MakeDelegate( this,
                    &Store<T,BY_AGE,BY_COHORT,BY_SPECIES,BY_GENOTYPE,BY_DRUG>::writeM), m );
        }
    }
    
    // Write stored values to stream
    void write( ostream& stream, size_t survey ){
        // use a (tree) map to sort by external measure
        map<int,size_t> measuresOrdered;
        for( size_t m = 0; m < outMeasures.size(); ++m ){
            measuresOrdered[outMeasures[m]] = m;
        }
        typedef pair<int,size_t> P;
        foreach( P mp, measuresOrdered ){
            writeM( stream, survey, mp.first, mp.second );
        }
    }
    
    // Checkpointing
    void checkpoint( ostream& stream ){
        reports.size() & stream;
        BOOST_FOREACH (T& y, reports) {
            y & stream;
        }
        // mIndices and outMeasures are constant after initialisation
    }
    void checkpoint( istream& stream ){
        size_t l;
        l & stream;
        if( l != size() ){
            throw util::checkpoint_error( "mon::reports: invalid list size" );
        }
        reports.resize (l);
        BOOST_FOREACH (T& y, reports) {
            y & stream;
        }
        // mIndices and outMeasures are constant after initialisation
    }
};

//NOTE: there may be more options than necessary. Optionally, A without C and
// C without A could be removed, and all outputs could be made doubles.
// Stores by integer value (no outputs include species or genotype):
Store<int, false, false, false, false, false> storeI;
Store<int, true, false, false, false, false> storeAI;
Store<int, false, true, false, false, false> storeCI;
Store<int, true, true, false, false, false> storeACI;
Store<int, false, false, false, true, false> storeGI;
Store<int, true, false, false, true, false> storeAGI;
Store<int, false, true, false, true, false> storeCGI;
Store<int, true, true, false, true, false> storeACGI;
Store<int, false, false, false, false, true> storePI;
Store<int, true, false, false, false, true> storeAPI;
Store<int, false, true, false, false, true> storeCPI;
Store<int, true, true, false, false, true> storeACPI;
// Stores by double value (note that by species reports never include age or cohort):
Store<double, false, false, false, false, false> storeF;
Store<double, true, false, false, false, false> storeAF;
Store<double, false, true, false, false, false> storeCF;
Store<double, true, true, false, false, false> storeACF;
Store<double, false, false, false, true, false> storeGF;
Store<double, true, false, false, true, false> storeAGF;
Store<double, false, true, false, true, false> storeCGF;
Store<double, true, true, false, true, false> storeACGF;
Store<double, false, false, false, false, true> storePF;
Store<double, true, false, false, false, true> storeAPF;
Store<double, false, true, false, false, true> storeCPF;
Store<double, true, true, false, false, true> storeACPF;
Store<double, false, false, true, false, false> storeSF;
Store<double, false, false, true, true, false> storeSGF;

int reportIMR = -1; // special output for fitting

void internal::initReporting( const scnXml::Scenario& scenario ){
    defineOutMeasures();
    const scnXml::MonitoringOptions& optsElt = scenario.getMonitoring().getSurveyOptions();
    set<int> outIds;    // all measure numbers used in output
    list<OutMeasure> enabledOutMeasures;
    foreach( const scnXml::MonitoringOption& optElt, optsElt.getOption() ){
        if( optElt.getValue() == false ) continue;      // option is disabled
        NamedMeasureMapT::const_iterator it =
            namedOutMeasures.find( optElt.getName() );
        if( it == namedOutMeasures.end() ){
            throw util::xml_scenario_error( (boost::format("unrecognised "
                "survey option: %1%") %optElt.getName()).str() );
        }
        OutMeasure om = it->second;     // copy; we may modify below
        if( om.m == M_NUM ){
            throw util::xml_scenario_error( (boost::format("obsolete "
                "survey option: %1%") %optElt.getName()).str() );
        }
        if( om.m == M_ALL_CAUSE_IMR ){
            if( om.isDouble && !om.byAge && !om.byCohort && !om.bySpecies ){
                reportIMR = om.outId;
            }else{
                throw util::xml_scenario_error( "measure allCauseIMR does not "
                    "support any categorisation" );
            }
        }
        if( optElt.getByAge().present() ){
            if( om.byAge ){
                om.byAge = optElt.getByAge().get();     // disable or keep
            }else if( optElt.getByAge().get() ){
                throw util::xml_scenario_error( (boost::format("measure %1% "
                    "does not support categorisation by age group")
                    %optElt.getName()).str() );
            }
        }
        if( optElt.getByCohort().present() ){
            if( om.byCohort ){
                om.byCohort = optElt.getByCohort().get();   // disable or keep
            }else if( optElt.getByCohort().get() ){
                throw util::xml_scenario_error( (boost::format("measure %1% "
                    "does not support categorisation by cohort")
                    %optElt.getName()).str() );
            }
        }
        if( optElt.getBySpecies().present() ){
            if( om.bySpecies ){
                om.bySpecies = optElt.getBySpecies().get(); // disable or keep
            }else if( optElt.getBySpecies().get() ){
                throw util::xml_scenario_error( (boost::format("measure %1% "
                    "does not support categorisation by age group")
                    %optElt.getName()).str() );
            }
        }
        if( optElt.getByGenotype().present() ){
            if( om.byGenotype ){
                om.byGenotype = optElt.getByGenotype().get();   // disable or keep
            }else if( optElt.getByGenotype().get() ){
                throw util::xml_scenario_error( (boost::format("measure %1% "
                    "does not support categorisation by age group")
                    %optElt.getName()).str() );
            }
        }
        if( optElt.getByDrugType().present() ){
            if( om.byDrug ){
                om.byDrug = optElt.getByDrugType().get();   // disable or keep
            }else if( optElt.getByDrugType().get() ){
                throw util::xml_scenario_error( (boost::format("measure %1% "
                    "does not support categorisation by drug type")
                    %optElt.getName()).str() );
            }
        }
        if( optElt.getOutputNumber().present() ) om.outId = optElt.getOutputNumber().get();
        if( outIds.count(om.outId) ){
            throw util::xml_scenario_error( (boost::format("monitoring output "
                "number %1% used more than once") %om.outId).str() );
        }
        outIds.insert( om.outId );
        enabledOutMeasures.push_back( om );
    }
    
    size_t nSpecies = scenario.getEntomology().getVector().present() ?
        scenario.getEntomology().getVector().get().getAnopheles().size() : 1;
    size_t nDrugs = scenario.getPharmacology().present() ?
        scenario.getPharmacology().get().getDrugs().getDrug().size() : 1;
    
    storeI.init( enabledOutMeasures, nSpecies, nDrugs );
    storeAI.init( enabledOutMeasures, nSpecies, nDrugs );
    storeCI.init( enabledOutMeasures, nSpecies, nDrugs );
    storeACI.init( enabledOutMeasures, nSpecies, nDrugs );
    storeGI.init( enabledOutMeasures, nSpecies, nDrugs );
    storeAGI.init( enabledOutMeasures, nSpecies, nDrugs );
    storeCGI.init( enabledOutMeasures, nSpecies, nDrugs );
    storeACGI.init( enabledOutMeasures, nSpecies, nDrugs );
    storePI.init( enabledOutMeasures, nSpecies, nDrugs );
    storeAPI.init( enabledOutMeasures, nSpecies, nDrugs );
    storeCPI.init( enabledOutMeasures, nSpecies, nDrugs );
    storeACPI.init( enabledOutMeasures, nSpecies, nDrugs );
    
    storeF.init( enabledOutMeasures, nSpecies, nDrugs );
    storeAF.init( enabledOutMeasures, nSpecies, nDrugs );
    storeCF.init( enabledOutMeasures, nSpecies, nDrugs );
    storeACF.init( enabledOutMeasures, nSpecies, nDrugs );
    storeGF.init( enabledOutMeasures, nSpecies, nDrugs );
    storeAGF.init( enabledOutMeasures, nSpecies, nDrugs );
    storeCGF.init( enabledOutMeasures, nSpecies, nDrugs );
    storeACGF.init( enabledOutMeasures, nSpecies, nDrugs );
    storePF.init( enabledOutMeasures, nSpecies, nDrugs );
    storeAPF.init( enabledOutMeasures, nSpecies, nDrugs );
    storeCPF.init( enabledOutMeasures, nSpecies, nDrugs );
    storeACPF.init( enabledOutMeasures, nSpecies, nDrugs );
    storeSF.init( enabledOutMeasures, nSpecies, nDrugs );
    storeSGF.init( enabledOutMeasures, nSpecies, nDrugs );
}

void internal::write( ostream& stream ){
    // use a (tree) map to sort by external measure
    typedef pair<WriteDelegate,size_t> MPair;
    map<int,MPair> measuresOrdered;
    
    storeI.addMeasures( measuresOrdered );
    storeAI.addMeasures( measuresOrdered );
    storeCI.addMeasures( measuresOrdered );
    storeACI.addMeasures( measuresOrdered );
    storeGI.addMeasures( measuresOrdered );
    storeAGI.addMeasures( measuresOrdered );
    storeCGI.addMeasures( measuresOrdered );
    storeACGI.addMeasures( measuresOrdered );
    storePI.addMeasures( measuresOrdered );
    storeAPI.addMeasures( measuresOrdered );
    storeCPI.addMeasures( measuresOrdered );
    storeACPI.addMeasures( measuresOrdered );
    
    storeF.addMeasures( measuresOrdered );
    storeAF.addMeasures( measuresOrdered );
    storeCF.addMeasures( measuresOrdered );
    storeACF.addMeasures( measuresOrdered );
    storeGF.addMeasures( measuresOrdered );
    storeAGF.addMeasures( measuresOrdered );
    storeCGF.addMeasures( measuresOrdered );
    storeACGF.addMeasures( measuresOrdered );
    storePF.addMeasures( measuresOrdered );
    storeAPF.addMeasures( measuresOrdered );
    storeCPF.addMeasures( measuresOrdered );
    storeACPF.addMeasures( measuresOrdered );
    storeSF.addMeasures( measuresOrdered );
    storeSGF.addMeasures( measuresOrdered );
    
    typedef pair<int,MPair> PP;
    for( size_t survey = 0; survey < impl::nSurveys; ++survey ){
        foreach( PP pp, measuresOrdered ){
            pp.second.first( stream, survey, pp.first, pp.second.second );
        }
    }
    if( reportIMR >= 0 ){
        // Infant mortality rate is a single number, therefore treated specially.
        // It is calculated across the entire intervention period and used in
        // model fitting.
        stream << 1 << "\t" << 1 << "\t" << reportIMR
            << "\t" << Clinical::infantAllCauseMort() << lineEnd;
    }
}

// Report functions: each reports to all usable stores (i.e. correct data type
// and where parameters don't have to be fabricated).
void reportMI( Measure measure, int val ){
    storeI.report( val, measure, impl::currentSurvey, 0, 0, 0, 0, 0 );
}
void reportMHI( Measure measure, const Host::Human& human, int val ){
    const size_t survey = impl::currentSurvey;
    const size_t ageIndex = human.monAgeGroup().i();
    storeI.report( val, measure, survey, 0, 0, 0, 0, 0 );
    storeAI.report( val, measure, survey, ageIndex, 0, 0, 0, 0 );
    storeCI.report( val, measure, survey, 0, human.cohortSet(), 0, 0, 0 );
    storeACI.report( val, measure, survey, ageIndex, human.cohortSet(), 0, 0, 0 );
}
void reportMSACI( Measure measure, size_t survey,
                  AgeGroup ageGroup, uint32_t cohortSet, int val )
{
    storeI.report( val, measure, survey, 0 ,0 ,0, 0, 0 );
    storeAI.report( val, measure, survey, ageGroup.i(), 0, 0, 0, 0 );
    storeCI.report( val, measure, survey, 0, cohortSet, 0, 0, 0 );
    storeACI.report( val, measure, survey, ageGroup.i(), cohortSet, 0, 0, 0 );
}
void reportMHGI( Measure measure, const Host::Human& human, size_t genotype,
                 int val )
{
    const size_t survey = impl::currentSurvey;
    const size_t ageIndex = human.monAgeGroup().i();
    storeI.report( val, measure, survey, 0, 0, 0, 0, 0 );
    storeAI.report( val, measure, survey, ageIndex, 0, 0, 0, 0 );
    storeCI.report( val, measure, survey, 0, human.cohortSet(), 0, 0, 0 );
    storeACI.report( val, measure, survey, ageIndex, human.cohortSet(), 0, 0, 0 );
    storeGI.report( val, measure, survey, 0, 0, 0, genotype, 0 );
    storeAGI.report( val, measure, survey, ageIndex, 0, 0, genotype, 0 );
    storeCGI.report( val, measure, survey, 0, human.cohortSet(), 0, genotype, 0 );
    storeACGI.report( val, measure, survey, ageIndex, human.cohortSet(), 0, genotype, 0 );
}
void reportMHPI( Measure measure, const Host::Human& human, size_t drugIndex,
                int val )
{
    const size_t survey = impl::currentSurvey;
    const size_t ageIndex = human.monAgeGroup().i();
    storeI.report( val, measure, survey, 0, 0, 0, 0, 0 );
    storeAI.report( val, measure, survey, ageIndex, 0, 0, 0, 0 );
    storeCI.report( val, measure, survey, 0, human.cohortSet(), 0, 0, 0 );
    storeACI.report( val, measure, survey, ageIndex, human.cohortSet(), 0, 0, 0 );
    storePI.report( val, measure, survey, 0, 0, 0, 0, drugIndex );
    storeAPI.report( val, measure, survey, ageIndex, 0, 0, 0, drugIndex );
    storeCPI.report( val, measure, survey, 0, human.cohortSet(), 0, 0, drugIndex );
    storeACPI.report( val, measure, survey, ageIndex, human.cohortSet(), 0, 0, drugIndex );
}
// Deployment reporting uses a different function to handle the method
// (mostly to make other types of report faster).
void reportMHD( Measure measure, const Host::Human& human,
                Deploy::Method method )
{
    const int val = 1;  // always report 1 deployment
    const size_t survey = impl::currentSurvey;
    size_t ageIndex = human.monAgeGroup().i();
    storeI.deploy( val, measure, survey, 0, 0, method );
    storeAI.deploy( val, measure, survey, ageIndex, 0, method );
    storeCI.deploy( val, measure, survey, 0, human.cohortSet(), method );
    storeACI.deploy( val, measure, survey, ageIndex, human.cohortSet(), method );
    // This is for nTreatDeployments:
    measure = MHD_ALL_DEPLOYS;
    storeI.deploy( val, measure, survey, 0, 0, method );
    storeAI.deploy( val, measure, survey, ageIndex, 0, method );
    storeCI.deploy( val, measure, survey, 0, human.cohortSet(), method );
    storeACI.deploy( val, measure, survey, ageIndex, human.cohortSet(), method );
}

void reportMF( Measure measure, double val ){
    storeF.report( val, measure, impl::currentSurvey, 0, 0, 0, 0, 0 );
}
void reportMHF( Measure measure, const Host::Human& human, double val ){
    const size_t survey = impl::currentSurvey;
    const size_t ageIndex = human.monAgeGroup().i();
    storeF.report( val, measure, survey, 0, 0, 0, 0, 0 );
    storeAF.report( val, measure, survey, ageIndex, 0, 0, 0, 0 );
    storeCF.report( val, measure, survey, 0, human.cohortSet(), 0, 0, 0 );
    storeACF.report( val, measure, survey, ageIndex, human.cohortSet(), 0, 0, 0 );
}
void reportMACGF( Measure measure, size_t ageIndex, uint32_t cohortSet,
                  size_t genotype, double val )
{
    const size_t survey = impl::currentSurvey;
    storeF.report( val, measure, survey, 0, 0, 0, 0, 0 );
    storeAF.report( val, measure, survey, ageIndex, 0, 0, 0, 0 );
    storeCF.report( val, measure, survey, 0, cohortSet, 0, 0, 0 );
    storeACF.report( val, measure, survey, ageIndex, cohortSet, 0, 0, 0 );
    storeGF.report( val, measure, survey, 0, 0, 0, genotype, 0 );
    storeAGF.report( val, measure, survey, ageIndex, 0, 0, genotype, 0 );
    storeCGF.report( val, measure, survey, 0, cohortSet, 0, genotype, 0 );
    storeACGF.report( val, measure, survey, ageIndex, cohortSet, 0, genotype, 0 );
}
void reportMHPF( Measure measure, const Host::Human& human, size_t drug, double val ){
    const size_t survey = impl::currentSurvey;
    const size_t ageIndex = human.monAgeGroup().i();
    storeF.report( val, measure, survey, 0, 0, 0, 0, 0 );
    storeAF.report( val, measure, survey, ageIndex, 0, 0, 0, 0 );
    storeCF.report( val, measure, survey, 0, human.cohortSet(), 0, 0, 0 );
    storeACF.report( val, measure, survey, ageIndex, human.cohortSet(), 0, 0, 0 );
    storePF.report( val, measure, survey, 0, 0, 0, 0, drug );
    storeAPF.report( val, measure, survey, ageIndex, 0, 0, 0, drug );
    storeCPF.report( val, measure, survey, 0, human.cohortSet(), 0, 0, drug );
    storeACPF.report( val, measure, survey, ageIndex, human.cohortSet(), 0, 0, drug );
}
void reportMHGF( Measure measure, const Host::Human& human, size_t genotype,
                 double val )
{
    reportMACGF( measure, human.monAgeGroup().i(), human.cohortSet(),
                 genotype, val );
}
void reportMSF( Measure measure, size_t species, double val ){
    const size_t survey = impl::currentSurvey;
    storeF.report( val, measure, survey, 0, 0, 0, 0, 0 );
    storeSF.report( val, measure, survey, 0, 0, species, 0, 0 );
}
void reportMSGF( Measure measure, size_t species, size_t genotype, double val ){
    const size_t survey = impl::currentSurvey;
    storeF.report( val, measure, survey, 0, 0, 0, 0, 0 );
    storeGF.report( val, measure, survey, 0, 0, 0, genotype, 0 );
    storeSF.report( val, measure, survey, 0, 0, species, 0, 0 );
    storeSGF.report( val, measure, survey, 0, 0, species, genotype, 0 );
}

bool isUsedM( Measure measure ){
    return storeI.isUsed(measure) ||
        storeAI.isUsed(measure) ||
        storeCI.isUsed(measure) ||
        storeACI.isUsed(measure) ||
        storeGI.isUsed(measure) ||
        storeAGI.isUsed(measure) ||
        storeCGI.isUsed(measure) ||
        storeACGI.isUsed(measure) ||
        storePI.isUsed(measure) ||
        storeAPI.isUsed(measure) ||
        storeCPI.isUsed(measure) ||
        storeACPI.isUsed(measure) ||
        storeF.isUsed(measure) ||
        storeAF.isUsed(measure) ||
        storeCF.isUsed(measure) ||
        storeACF.isUsed(measure) ||
        storeGF.isUsed(measure) ||
        storeAGF.isUsed(measure) ||
        storeCGF.isUsed(measure) ||
        storeACGF.isUsed(measure) ||
        storePF.isUsed(measure) ||
        storeAPF.isUsed(measure) ||
        storeCPF.isUsed(measure) ||
        storeACPF.isUsed(measure) ||
        storeSF.isUsed(measure) ||
        storeSGF.isUsed(measure);
}

void checkpoint( ostream& stream ){
    impl::currentSurvey & stream;
    impl::nextSurveyTime & stream;
    
    storeI.checkpoint(stream);
    storeAI.checkpoint(stream);
    storeCI.checkpoint(stream);
    storeACI.checkpoint(stream);
    storeGI.checkpoint(stream);
    storeAGI.checkpoint(stream);
    storeCGI.checkpoint(stream);
    storeACGI.checkpoint(stream);
    storePI.checkpoint(stream);
    storeAPI.checkpoint(stream);
    storeCPI.checkpoint(stream);
    storeACPI.checkpoint(stream);
    
    storeF.checkpoint(stream);
    storeAF.checkpoint(stream);
    storeCF.checkpoint(stream);
    storeACF.checkpoint(stream);
    storeGF.checkpoint(stream);
    storeAGF.checkpoint(stream);
    storeCGF.checkpoint(stream);
    storeACGF.checkpoint(stream);
    storePF.checkpoint(stream);
    storeAPF.checkpoint(stream);
    storeCPF.checkpoint(stream);
    storeACPF.checkpoint(stream);
    storeSF.checkpoint(stream);
    storeSGF.checkpoint(stream);
}
void checkpoint( istream& stream ){
    impl::currentSurvey & stream;
    impl::nextSurveyTime & stream;
    
    storeI.checkpoint(stream);
    storeAI.checkpoint(stream);
    storeCI.checkpoint(stream);
    storeACI.checkpoint(stream);
    storeGI.checkpoint(stream);
    storeAGI.checkpoint(stream);
    storeCGI.checkpoint(stream);
    storeACGI.checkpoint(stream);
    storePI.checkpoint(stream);
    storeAPI.checkpoint(stream);
    storeCPI.checkpoint(stream);
    storeACPI.checkpoint(stream);
    
    storeF.checkpoint(stream);
    storeAF.checkpoint(stream);
    storeCF.checkpoint(stream);
    storeACF.checkpoint(stream);
    storeGF.checkpoint(stream);
    storeAGF.checkpoint(stream);
    storeCGF.checkpoint(stream);
    storeACGF.checkpoint(stream);
    storePF.checkpoint(stream);
    storeAPF.checkpoint(stream);
    storeCPF.checkpoint(stream);
    storeACPF.checkpoint(stream);
    storeSF.checkpoint(stream);
    storeSGF.checkpoint(stream);
}

}
}
