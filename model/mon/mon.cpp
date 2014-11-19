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

#include "mon/reporting.h"
#include "Monitoring/Survey.h"
#include "Monitoring/Surveys.h"
#include "Host/Human.h"
#define H_OM_mon_cpp
#include "mon/OutputMeasures.hpp"
#include "schema/monitoring.h"

#include <FastDelegate.h>
#include <iostream>
#include <boost/format.hpp>

namespace OM {
namespace mon {
    using namespace fastdelegate;

// Constants or defined during init:
size_t nSurveys = 0;
// Accumulators, variables:
size_t currentSurvey = NOT_USED;

#ifndef NDEBUG
const size_t NOT_ACCEPTED = boost::integer_traits<size_t>::const_max - 1;
#endif

typedef FastDelegate4<ostream&,size_t,int,size_t> WriteDelegate;

// Store by human age and cohort
template<typename T, bool BY_AGE, bool BY_COHORT, bool BY_SPECIES>
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
    
    size_t nAgeGroups, nCohortSets, nSpecies;
    // These are the stored reports
    // Array has four dimensions; see index()
    vector<T> reports;
    
    // get size of reports
    inline size_t size(){ return outMeasures.size() * nSurveys *
        nAgeGroups * nCohortSets * nSpecies; }
    // get an index in reports
    inline size_t index( size_t m, size_t s, size_t a, size_t c, size_t sp ){
        return sp + nSpecies * (c + nCohortSets * (a + nAgeGroups * (s + nSurveys * m)));
    }
    
    inline void add(size_t mIndex, size_t survey, size_t ageIndex,
                    uint32_t cohortSet, size_t species, T val)
    {
#ifndef NDEBUG
        if( mIndex >= outMeasures.size() ||
            survey >= nSurveys ||
            ageIndex >= nAgeGroups ||
            cohortSet >= nCohortSets ||
            species >= nSpecies
        ){
            cout << "Index out of bounds for survey:\t" << survey << " of " << nSurveys
                << "\nmeasure number\t" << mIndex << " of " << outMeasures.size()
                << "\nage group\t" << ageIndex << " of " << nAgeGroups
                << "\ncohort set\t" << cohortSet << " of " << nCohortSets
                << "\nspecies\t" << species << " of " << nSpecies
                << endl;
        }
#endif
        reports[index(mIndex,survey,ageIndex,cohortSet,species)] += val;
    }
    
    void writeM( ostream& stream, size_t survey, int outMeasure, size_t inMeasure ){
        if( BY_SPECIES ){
            assert( !BY_AGE && !BY_COHORT );  // output col2 conflicts
            for( size_t species = 0; species < nSpecies; ++species ){
                //TODO: should we make this backwards-compatible (output
                // name, not a number)?
                const int col2 = species + 1;
                T value = reports[index(inMeasure,survey,0,0,species)];
                stream << (survey+1) << '\t' << col2 << '\t' << outMeasure
                    << '\t' << value << Monitoring::lineEnd;
            }
        }else{
            // Backwards compatibility: first age group starts at 1, unless
            // there isn't an age group:
            const int ageGroupAdd = BY_AGE ? 1 : 0;
            for( size_t cohortSet = 0; cohortSet < nCohortSets; ++cohortSet ){
            for( size_t ageGroup = 0; ageGroup < nAgeGroups; ++ageGroup ){
                // Yeah, >999 age groups clashes with cohort sets, but unlikely a real issue
                const int col2 = ageGroup + ageGroupAdd +
                    1000 * Monitoring::Surveys.cohortSetOutputId( cohortSet );
                T value = reports[index(inMeasure,survey,ageGroup,cohortSet,0)];
                stream << (survey+1) << '\t' << col2 << '\t' << outMeasure
                    << '\t' << value << Monitoring::lineEnd;
            } }
        }
    }
    
public:
    // Set up ready to accept reports.
    void init( const list<OutMeasure>& required, size_t nAG, size_t nCS, size_t nSp ){
        nAgeGroups = BY_AGE ? nAG : 1;
        nCohortSets = BY_COHORT ? nCS : 1;
        nSpecies = (BY_SPECIES && nSp > 0) ? nSp : 1;
        mIndices.assign( M_NUM, NOT_USED );
        // outMeasures.size() is the number of measures we store here
        outMeasures.assign( 0, 0 );
        
        for( list<OutMeasure>::const_iterator it = required.begin();
            it != required.end(); ++it )
        {
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
                it->bySpecies != BY_SPECIES ) continue;
            
            if( mIndices[it->m] != NOT_USED ||
                (it->method == Deploy::NA && deployIndices.count(it->m) > 0) )
            {
                //FIXME: use a name not a number to describe the measure
                throw util::xml_scenario_error( (boost::format("multiple use "
                    "of output measure %1% by age and cohort") %it->m).str() );
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
    void report( T val, Measure measure, size_t survey, size_t ageIndex = 0,
                 uint32_t cohortSet = 0, size_t species = 0 )
    {
        if( survey == NOT_USED ) return; // pre-main-sim & unit tests we ignore all reports
        assert( mIndices[measure] != NOT_ACCEPTED );
        //TODO: is this optimisable?
        assert( BY_AGE || (nAgeGroups == 1 && ageIndex == 0) );
        assert( BY_COHORT || (nCohortSets == 1 && cohortSet == 0) );
        assert( BY_SPECIES || (nSpecies == 1 && species == 0) );
        if( mIndices[measure] == NOT_USED ){    // measure not used by this store
            assert( deployIndices.count(measure) == 0 );
            return;
        }
        // last category is for humans too old for reporting groups:
        if( ageIndex == nAgeGroups ) return;
        add( mIndices[measure], survey, ageIndex, cohortSet, species, val );
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
        assert( mIndices[measure] == NOT_USED );
        assert( nSpecies == 1 );        // this function does not take species
        
        typedef DeployInd_t::const_iterator const_it_t;
        pair<const_it_t,const_it_t> range = deployIndices.equal_range( measure );
        for( const_it_t it = range.first; it != range.second; ++it ){
            uint8_t mask = it->second.first;
            if( (mask & method) == 0 ) continue;
            add( it->second.second, survey, ageIndex, cohortSet, 0, val );
        }
    }
    
    // Order self in a list of outputs
    void addMeasures( map<int,pair<WriteDelegate,size_t> >& mOrdered ){
        for( size_t m = 0; m < outMeasures.size(); ++m ){
            mOrdered[outMeasures[m]] = make_pair( MakeDelegate( this,
                    &Store<T,BY_AGE,BY_COHORT,BY_SPECIES>::writeM), m );
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

Store<int, false, false, false> storeI;
Store<double, false, false, false> storeF;
Store<int, true, true, false> storeHACI;
Store<double, true, true, false> storeHACF;
Store<double, false, false, true> storeSF;

void initialise( size_t numSurveys, size_t nCohortSets, size_t nSpecies,
                 const scnXml::Monitoring& monElt )
{
    nSurveys = numSurveys;
    // Last group is those too old to be reported:
    size_t nAgeGroups = Monitoring::AgeGroup::getNumGroups() - 1;
    currentSurvey = NOT_USED;
    
    defineOutMeasures();
    const scnXml::OptionSet& optsElt = monElt.getSurveyOptions();
    list<OutMeasure> enabledOutMeasures;
    foreach( const scnXml::Option& optElt, optsElt.getOption() ){
        if( optElt.getValue() == false ) continue;      // option is disabled
        NamedMeasureMapT::const_iterator it = namedOutMeasures.find( optElt.getName() );
        if( it == namedOutMeasures.end() ) continue;    //TODO: eventually throw an xml_scenario_error
        enabledOutMeasures.push_back( it->second );
    }
    storeI.init( enabledOutMeasures, nAgeGroups, nCohortSets, nSpecies );
    storeF.init( enabledOutMeasures, nAgeGroups, nCohortSets, nSpecies );
    storeSF.init( enabledOutMeasures, nAgeGroups, nCohortSets, nSpecies );
    storeHACI.init( enabledOutMeasures, nAgeGroups, nCohortSets, nSpecies );
    storeHACF.init( enabledOutMeasures, nAgeGroups, nCohortSets, nSpecies );
}

void initMainSim(){
    currentSurvey = 0;
}
void concludeSurvey(){
    currentSurvey += 1;
    // After the last survey has completed:
    if( currentSurvey >= nSurveys ) currentSurvey = NOT_USED;
}

void write( ostream& stream ){
    // use a (tree) map to sort by external measure
    typedef pair<WriteDelegate,size_t> MPair;
    map<int,MPair> measuresOrdered;
    storeI.addMeasures( measuresOrdered );
    storeF.addMeasures( measuresOrdered );
    storeSF.addMeasures( measuresOrdered );
    storeHACI.addMeasures( measuresOrdered );
    storeHACF.addMeasures( measuresOrdered );
    
    typedef pair<int,MPair> PP;
    for( size_t survey = 0; survey < nSurveys; ++survey ){
        foreach( PP pp, measuresOrdered ){
            pp.second.first( stream, survey, pp.first, pp.second.second );
        }
    }
}

// Report functions: each reports to all usable stores (i.e. correct data type
// and where parameters don't have to be fabricated).
void reportMI( Measure measure, int val ){
    storeI.report( val, measure, currentSurvey );
}
void reportMF( Measure measure, double val ){
    storeF.report( val, measure, currentSurvey );
}
void reportMHI( Measure measure, const Host::Human& human, int val ){
    storeI.report( val, measure, currentSurvey );
    size_t ageIndex = human.getMonitoringAgeGroup().i();
    storeHACI.report( val, measure, currentSurvey, ageIndex, human.cohortSet() );
}
void reportMHF( Measure measure, const Host::Human& human, double val ){
    storeF.report( val, measure, currentSurvey );
    size_t ageIndex = human.getMonitoringAgeGroup().i();
    storeHACF.report( val, measure, currentSurvey, ageIndex, human.cohortSet() );
}
void reportMSACI( Measure measure, size_t survey,
                  Monitoring::AgeGroup ageGroup, uint32_t cohortSet, int val )
{
    storeI.report( val, measure, currentSurvey );
    storeHACI.report( val, measure, survey, ageGroup.i(), cohortSet );
}
void reportMACF( Measure measure, size_t ageIndex, uint32_t cohortSet, double val )
{
    storeF.report( val, measure, currentSurvey );
    storeHACF.report( val, measure, currentSurvey, ageIndex, cohortSet );
}
void reportMSF( Measure measure, size_t species, double val ){
    storeF.report( val, measure, currentSurvey );
    storeSF.report( val, measure, currentSurvey, 0, 0, species );
}
// Deployment reporting uses a different function to handle the method
// (mostly to make other types of report faster).
void reportMHD( Measure measure, const Host::Human& human,
                Deploy::Method method )
{
    const int val = 1;  // always report 1 deployment
    storeI.deploy( val, measure, currentSurvey, 0, 0, method );
    size_t ageIndex = human.getMonitoringAgeGroup().i();
    storeHACI.deploy( val, measure, currentSurvey, ageIndex, human.cohortSet(), method );
    // This is for nTreatDeployments:
    storeHACI.deploy( val, MHD_ALL_DEPLOYS, currentSurvey, ageIndex,
                      human.cohortSet(), method );
}

void checkpoint( ostream& stream ){
    currentSurvey & stream;
    storeI.checkpoint(stream);
    storeF.checkpoint(stream);
    storeSF.checkpoint(stream);
    storeHACI.checkpoint(stream);
    storeHACF.checkpoint(stream);
}
void checkpoint( istream& stream ){
    currentSurvey & stream;
    storeI.checkpoint(stream);
    storeF.checkpoint(stream);
    storeSF.checkpoint(stream);
    storeHACI.checkpoint(stream);
    storeHACF.checkpoint(stream);
}

}
}
