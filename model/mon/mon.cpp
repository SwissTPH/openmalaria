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

#include <iostream>
#include <boost/format.hpp>

namespace OM {
namespace mon {

// Constants or defined during init:
size_t nSurveys, nAgeGroups, nCohortSets;
// Accumulators, variables:
size_t currentSurvey;

#ifndef NDEBUG
const size_t NOT_ACCEPTED = boost::integer_traits<size_t>::const_max - 1;
#endif

// Store by human age and cohort
template<typename T, bool BY_AGE, bool BY_COHORT>
struct Store{
    // These are the stored reports
    // Array has four dimensions; see index()
    vector<T> reports;
    // get an index in reports
    inline size_t index( size_t m, size_t s, size_t a, size_t c ){
        size_t i = s + nSurveys * m;
        if( BY_AGE ) i = a + nAgeGroups * i;
        if( BY_COHORT ) i = c + nCohortSets * i;
        return i;
    }
    
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
    
    // Set up ready to accept reports.
    void init( const list<OutMeasure>& required ){
        mIndices.assign( M_NUM, NOT_USED );
        // outMeasures.size() is the number of measures we store here
        outMeasures.assign( 0, 0 );
        
        for( list<OutMeasure>::const_iterator it = required.begin(); it != required.end(); ++it ){
            if( it->isDouble != (typeid(T) == typeid(double) ) ){
#ifndef NDEBUG
                // Debug mode: this should prevent silly errors where the type reported does not match the type defined for some output
                mIndices[it->m] = NOT_ACCEPTED;
#endif
                continue;
            }
            if( it->byAge != BY_AGE || it->byCohort != BY_COHORT ) continue;
            
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
        
        reports.assign(outMeasures.size() * nSurveys * nAgeGroups * nCohortSets, 0);
    }
    
    // Take a reported value and either store it or forget it.
    // If some ageIndex or cohortSet are not applicable, use 0.
    void report( Measure measure, size_t survey, size_t ageIndex, uint32_t cohortSet, T val ){
        assert( mIndices[measure] != NOT_ACCEPTED );
        if( survey == NOT_USED ) return; // pre-main-sim we ignore all reports
        if( mIndices[measure] == NOT_USED ){    // measure not used by this store
            assert( deployIndices.count(measure) == 0 );
            return;
        }
        if( ageIndex == nAgeGroups ) return;    // last category is for humans too old for reporting groups
        add( mIndices[measure], survey, ageIndex, cohortSet, val );
    }
    
    // Take a deployment report and potentially store it in one or more places
    // If some ageIndex or cohortSet are not applicable, use 0.
    void deploy( Measure measure, size_t survey, size_t ageIndex, uint32_t cohortSet, Deploy::Method method, T val ){
        assert( method == Deploy::NA || method == Deploy::TIMED || method == Deploy::CTS || method == Deploy::TREAT );
        assert( mIndices[measure] == NOT_USED );
        if( survey == NOT_USED ) return; // pre-main-sim we ignore all reports
        if( ageIndex == nAgeGroups ) return;    // last category is for humans too old for reporting groups
        
        typedef DeployInd_t::const_iterator const_it_t;
        pair<const_it_t,const_it_t> range = deployIndices.equal_range( measure );
        for( const_it_t it = range.first; it != range.second; ++it ){
            uint8_t mask = it->second.first;
            if( (mask & method) == 0 ) continue;
            add( it->second.second, survey, ageIndex, cohortSet, val );
        }
    }
    
    // Write stored values to stream
    void write( ostream& stream, int survey ){
        map<int,size_t> measuresOrdered;
        for( size_t m = 0; m < outMeasures.size(); ++m ){
            measuresOrdered[outMeasures[m]] = m;
        }
        typedef pair<int,size_t> P;
        size_t nCohorts = BY_COHORT ? nCohortSets : 1;
        size_t nAges = BY_AGE ? nAgeGroups : 1;
        const int addCol2 = BY_AGE ? 1 : 0;     // backwards compatibility: first age group starts at 0, unless there isn't an age group
        foreach( P mp, measuresOrdered ){
            for( size_t cohortSet = 0; cohortSet < nCohorts; ++cohortSet ){
                for( size_t ageGroup = 0; ageGroup < nAges; ++ageGroup ){
                    // Yeah, >999 age groups clashes with cohort sets, but unlikely a real issue
                    const int col2 = 1000 * Monitoring::Surveys.cohortSetOutputId( cohortSet ) + ageGroup + addCol2;
                    T value = reports[index(mp.second,survey-1,ageGroup,cohortSet)];
                    stream << survey << '\t' << col2 << '\t' << mp.first
                        << '\t' << value << Monitoring::lineEnd;
                }
            }
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
        if( l != outMeasures.size() * nSurveys * nAgeGroups * nCohortSets ){
            throw util::checkpoint_error( "mon::reports: invalid list size" );
        }
        reports.resize (l);
        BOOST_FOREACH (T& y, reports) {
            y & stream;
        }
        // mIndices and outMeasures are constant after initialisation
    }
    
private:
    inline void add(size_t mIndex, size_t survey, size_t ageIndex, uint32_t cohortSet, T val){
#ifndef NDEBUG
        if( mIndex >= outMeasures.size() ||
            survey >= nSurveys ||
            ageIndex >= nAgeGroups ||
            cohortSet >= nCohortSets
        ){
            cout << "Index out of bounds for survey:\t" << survey << " of " << nSurveys
                << "\nmeasure number\t" << mIndex << " of " << outMeasures.size()
                << "\nage group\t" << ageIndex << " of " << nAgeGroups
                << "\ncohort set\t" << cohortSet << " of " << nCohortSets
                << endl;
        }
#endif
        reports[index(mIndex,survey,ageIndex,cohortSet)] += val;
    }
};

Store<int, false, false> storeI;
Store<double, false, false> storeF;
Store<int, true, true> storeHACI;
Store<double, true, true> storeHACF;

void initialise( size_t nS, size_t nAG, size_t nCS, const scnXml::Monitoring& monElt ){
    nSurveys = nS - 1;  //NOTE: nS is actually 1 greater than required
    nAgeGroups = nAG - 1;       //NOTE: last group is those too old to be reported
    nCohortSets = nCS;
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
    storeI.init( enabledOutMeasures );
    storeF.init( enabledOutMeasures );
    storeHACI.init( enabledOutMeasures );
    storeHACF.init( enabledOutMeasures );
}

void initMainSim(){
    currentSurvey = 0;
}
void concludeSurvey(){
    currentSurvey += 1;
    // After the last survey has completed:
    if( currentSurvey >= nSurveys ) currentSurvey = NOT_USED;
}

void write( ostream& stream, int survey ){
    //TODO: re-order these by output measure
    storeHACI.write( stream, survey );
    storeHACF.write( stream, survey );
    storeI.write( stream, survey );
    storeF.write( stream, survey );
}

// Report functions: each reports to all relevant stores.
void reportMI( Measure measure, int val ){
    storeI.report( measure, currentSurvey, 0, 0, val );
}
void reportMF( Measure measure, double val ){
    storeF.report( measure, currentSurvey, 0, 0, val );
}
void reportMHI( Measure measure, const Host::Human& human, int val ){
    storeI.report( measure, currentSurvey, 0, 0, val );
    size_t ageIndex = human.getMonitoringAgeGroup().i();
    storeHACI.report( measure, currentSurvey, ageIndex, human.cohortSet(), val );
}
void reportMHF( Measure measure, const Host::Human& human, double val ){
    storeF.report( measure, currentSurvey, 0, 0, val );
    size_t ageIndex = human.getMonitoringAgeGroup().i();
    storeHACF.report( measure, currentSurvey, ageIndex, human.cohortSet(), val );
}
void reportMSACI( Measure measure, size_t survey, Monitoring::AgeGroup ageGroup, uint32_t cohortSet, int val ){
    storeI.report( measure, currentSurvey, 0, 0, val );
    storeHACI.report( measure, survey, ageGroup.i(), cohortSet, val );
}
// Deployment reporting uses a different function to handle the method
// (mostly to make other types of report faster).
void reportMHD( Measure measure, const Host::Human& human, Deploy::Method method ){
    const int val = 1;  // always report 1 deployment
    storeI.deploy( measure, currentSurvey, 0, 0, method, val );
    size_t ageIndex = human.getMonitoringAgeGroup().i();
    storeHACI.deploy( measure, currentSurvey, ageIndex, human.cohortSet(), method, val );
    // This is for nTreatDeployments:
    storeHACI.deploy( MHD_ALL_DEPLOYS, currentSurvey, ageIndex, human.cohortSet(), method, val );
}

void checkpoint( ostream& stream ){
    currentSurvey & stream;
    storeI.checkpoint(stream);
    storeF.checkpoint(stream);
    storeHACI.checkpoint(stream);
    storeHACF.checkpoint(stream);
}
void checkpoint( istream& stream ){
    currentSurvey & stream;
    storeI.checkpoint(stream);
    storeF.checkpoint(stream);
    storeHACI.checkpoint(stream);
    storeHACF.checkpoint(stream);
}

}
}
