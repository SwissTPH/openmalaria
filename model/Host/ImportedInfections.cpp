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

#include "Host/ImportedInfections.h"
#include "Host/Human.h"
#include "util/random.h"
#include "util/UnitParse.h"
#include "Population.h"

namespace OM { namespace Host {

void ImportedInfections::init( const scnXml::ImportedInfections& iiElt ){
    const scnXml::ImportedInfections::TimedType& tElt = iiElt.getTimed();
    try{
        //NOTE: if changing XSD, this should not have a default unit:
        period = UnitParse::readDuration( tElt.getPeriod(), UnitParse::STEPS );
        if( period < SimTime::zero() ){
            throw util::format_error( "cannot be negative" );
        }
    }catch( const util::format_error& e ){
        throw util::xml_scenario_error( string("interventions/importedInfections/timed/period: ").append(e.message()) );
    }
    rate.reserve( tElt.getRate().size() );
    try{
        for( auto it = tElt.getRate().begin(); it != tElt.getRate().end(); ++it ){
            SimDate date = UnitParse::readDate( it->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/ );
            // convert to per-time-step, per-person
            double rateVal = it->getValue() * sim::yearsPerStep() * (1.0 / 1000.0);
            rate.push_back( Rate( date - sim::startDate(), rateVal ) );
        }
        sort( rate.begin(), rate.end() );
        if( rate.size() > 0 ){
            if( period != SimTime::zero() && rate[0].time != SimTime::zero() ){
                throw util::xml_scenario_error( "must specify rate at time zero when period is not zero" );
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
        }
    }catch( const util::format_error& e ){
        throw util::xml_scenario_error( string("interventions/importedInfections/timed/time: ").append(e.message()) );
    }
}

void ImportedInfections::import( Population& population ){
    if( rate.size() == 0 ) return;      // no imported infections
    SimTime now = sim::intervTime();
    assert( now >= SimTime::zero() );
    if( period > SimTime::zero() ){
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
        for(Human& human : population){
            if(human.rng().bernoulli( rateNow )){
                human.addInfection();
            }
        }
    }
}

} }
