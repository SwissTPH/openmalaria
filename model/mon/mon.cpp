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
template<typename T>
struct StoreHAC{
    // These are the stored reports
    // Array has four dimensions; see index()
    vector<T> reports;
    // This maps from measures (MHR_HOSTS, etc) to an index in reports or NOT_USED
    vector<size_t> mIndices;
    // This maps from an index in reports to an output measure
    vector<int> outMeasures;
    
    // get an index in reports
    inline size_t index( size_t m, size_t s, size_t a, size_t c ){
        return c + nCohortSets * (a + nAgeGroups * (s + nSurveys * m));
    }
    
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
            if( !it->byAge || !it->byCohort ) continue;
            
            if( mIndices[it->m] != NOT_USED ){
                //FIXME: use a name not a number to describe the measure
                throw util::xml_scenario_error( (boost::format("multiple use "
                    "of output measure %1% by age and cohort") %it->m).str() );
            }
            
            mIndices[it->m] = outMeasures.size();       // length becomes our index
            outMeasures.push_back( it->outId ); // increment length
        }
        
        reports.resize(outMeasures.size() * nSurveys * nAgeGroups * nCohortSets);
    }
    
    // Take a value and either store it or forget it
    void accept( Measure measure, size_t survey, size_t ageIndex, uint32_t cohortSet, T val ){
        if( survey == NOT_USED ) return; // pre-main-sim we ignore all reports
        assert( mIndices[measure] != NOT_ACCEPTED );
        if( mIndices[measure] == NOT_USED ) return;     // measure not used by this store
        if( ageIndex == nAgeGroups ) return;    // last category is for humans too old for reporting groups
        size_t mIndex = mIndices[measure];
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
    
    // Write stored values to stream
    void write( ostream& stream, int survey ){
        for( size_t m = 0; m < outMeasures.size(); ++m ){
            int outMeasure = outMeasures[m];
            for( size_t cohortSet = 0; cohortSet < nCohortSets; ++cohortSet ){
                for( size_t ageGroup = 0; ageGroup < nAgeGroups; ++ageGroup ){
                    // Yeah, >999 age groups clashes with cohort sets, but unlikely a real issue
                    const int col2 = 1000 * Monitoring::Surveys.cohortSetOutputId( cohortSet ) + ageGroup + 1;
                    T value = reports[index(m,survey-1,ageGroup,cohortSet)];
                    stream << survey << '\t' << col2 << '\t' << outMeasure
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
};

StoreHAC<int> storeHACI;
StoreHAC<double> storeHACD;

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
    storeHACI.init( enabledOutMeasures );
    storeHACD.init( enabledOutMeasures );
}

void initMainSim(){
    currentSurvey = 0;
}
void concludeSurvey(){
    currentSurvey += 1;
    // After the last survey has completed:
    if( currentSurvey >= nSurveys ) currentSurvey = NOT_USED;
}

void writeMHI( ostream& stream, int survey ){
    storeHACI.write( stream, survey );
}
void writeMHD( ostream& stream, int survey ){
    storeHACD.write( stream, survey );
}

void reportMHI( Measure measure, const Host::Human& human, int val ){
    size_t ageIndex = human.getMonitoringAgeGroup().i();
    storeHACI.accept( measure, currentSurvey, ageIndex, human.cohortSet(), val );
}
void reportMHD( Measure measure, const Host::Human& human, double val ){
    size_t ageIndex = human.getMonitoringAgeGroup().i();
    storeHACD.accept( measure, currentSurvey, ageIndex, human.cohortSet(), val );
}
void reportMSACI( Measure measure, size_t survey, Monitoring::AgeGroup ageGroup, uint32_t cohortSet, int val ){
    storeHACI.accept( measure, survey, ageGroup.i(), cohortSet, val );
}

void checkpoint( ostream& stream ){
    currentSurvey & stream;
    storeHACI.checkpoint(stream);
    storeHACD.checkpoint(stream);
}
void checkpoint( istream& stream ){
    currentSurvey & stream;
    storeHACI.checkpoint(stream);
    storeHACD.checkpoint(stream);
}

}
}
