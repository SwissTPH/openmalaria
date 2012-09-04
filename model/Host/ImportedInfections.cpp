/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2011 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#include "Host/ImportedInfections.h"
#include "util/random.h"

namespace OM { namespace Host {
    
bool ImportedInfections::init( const scnXml::ImportedInfections& iiElt ){
    const scnXml::ImportedInfections::TimedType& tElt = iiElt.getTimed();
    if( tElt.getPeriod() < 0 ){
        throw util::xml_scenario_error( "interventions.importedInfections.timed.period cannot be negative" );
    }
    period = TimeStep( tElt.getPeriod() );
    rate.reserve( tElt.getRate().size() );
    for( scnXml::ImportedInfections::TimedType::RateSequence::const_iterator it = tElt.getRate().begin(); it != tElt.getRate().end(); ++it ){
        TimeStep time( it->getTime() );
        if( period != TimeStep(0) && time >= period ){
            throw util::xml_scenario_error( "interventions.importedInfections.timed: time cannot be greater than period when period is not zero" );
        }
        // convert to per-timestep, per-person
        double rateVal = it->getValue() / TimeStep::intervalsPerYear.asInt() / 1000.0;
        rate.push_back( Rate( time, rateVal ) );
    }
    sort( rate.begin(), rate.end() );
    if( rate.size() > 0 ){
        if( period != TimeStep(0) && rate[0].time != TimeStep( 0 ) ){
            throw util::xml_scenario_error( "interventions.importedInfections.timed: must specify rate at time zero when period is not zero" );
        }
        // remove useless repeated entries from list
        double lastRate = rate[0].value;
        for( size_t i = 1; i < rate.size(); ){
            if( rate[i].value == lastRate ){
                rate.erase( rate.begin() + i );
            }else{
                lastRate = rate[i].value;
                ++i;
            }
        }
        // first or second value non-zero
        return (rate[0].value != 0.0) || (rate.size() > 1 && rate[1].value != 0.0);
    }
    return false;       // no list
}

void ImportedInfections::import( Population& population ){
    assert( rate.size() > 0 );  // please don't call me otherwise!
    assert( TimeStep::interventionPeriod >= TimeStep(0) );
    TimeStep now = TimeStep::interventionPeriod;
    if( period > TimeStep(0) ){
        now = mod_nn(now, period);
    }
    if( rate[lastIndex].time > now ){
        lastIndex = 0;  // gone round in a loop: back to start of period
    }
    if( now > rate[lastIndex].time && rate.size() > lastIndex+1 ){
        if( rate[lastIndex+1].time <= now ){
            ++lastIndex;
            // won't be necessary to do this more than once, since import() is called every time step
        }
    }
    
    double rateNow = rate[lastIndex].value;
    if( rateNow > 0.0 ){
        for(Population::HumanIter it = population.getList().begin(); it!=population.getList().end(); ++it){
            if(util::random::bernoulli( rateNow )){
                it->addInfection();
            }
        }
    }
}

} }
