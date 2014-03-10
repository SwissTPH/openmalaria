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

#include "Clinical/Episode.h"
#include "Monitoring/Surveys.h"

namespace OM {
namespace Clinical {
using Monitoring::Surveys;
using Monitoring::Survey;

TimeStep Episode::healthSystemMemory( TimeStep::never );


void Episode::init( int hsMemory ) {
    healthSystemMemory = TimeStep( hsMemory );
}


Episode::~Episode ()
{
    report ();
}

void Episode::flush() {
    report();
    _time = TimeStep::never;
}


void Episode::update (bool inCohort, Monitoring::AgeGroup ageGroup, Episode::State newState)
{
    if (TimeStep::simulation > (_time + healthSystemMemory)) {
        report ();

        _time = TimeStep::simulation;
        _surveyPeriod = Surveys.getSurveyNumber( inCohort );
        _ageGroup = ageGroup;
        _state = newState;
    } else {
        _state = Episode::State (_state | newState);
    }
}

void Episode::report () {
    if (_time == TimeStep::never)        // Nothing to report
        return;

    // Reports malarial/non-malarial UC fever dependent on cause, not diagnosis.
    if (_state & Episode::MALARIA) {
        // Malarial fevers: report bout
        if (_state & Episode::COMPLICATED) {
            Surveys.at(_surveyPeriod)
                .addInt( Survey::MI_SEVERE_EPISODES, _ageGroup, 1) ;
        } else { // UC or UC2
            Surveys.at(_surveyPeriod)
                .addInt( Survey::MI_UNCOMPLICATED_EPISODES, _ageGroup, 1 );
        }

        // Report outcomes of malarial fevers
        if (_state & Episode::EVENT_IN_HOSPITAL) {
            if (_state & Episode::DIRECT_DEATH) {
                Surveys.at(_surveyPeriod)
                    .addInt( Survey::MI_DIRECT_DEATHS, _ageGroup, 1 )
                    .addInt( Survey::MI_HOSPITAL_DEATHS, _ageGroup, 1 );
                if (_state & Episode::EVENT_FIRST_DAY)
                    Surveys.at(_surveyPeriod)
                        .addInt( Survey::MI_FIRST_DAY_DEATHS, _ageGroup, 1 )
                        .addInt( Survey::MI_HOSPITAL_FIRST_DAY_DEATHS, _ageGroup, 1 );
            }
            else if (_state & Episode::SEQUELAE) {
                Surveys.at(_surveyPeriod)
                    .addInt( Survey::MI_SEQUELAE, _ageGroup, 1 )
                    .addInt( Survey::MI_HOSPITAL_SEQUELAE, _ageGroup, 1 );
            }
            else if (_state & Episode::RECOVERY)
                Surveys.at(_surveyPeriod)
                    .addInt( Survey::MI_HOSPITAL_RECOVERIES, _ageGroup, 1 );
        } else {
            if (_state & Episode::DIRECT_DEATH) {
                Surveys.at(_surveyPeriod)
                    .addInt( Survey::MI_DIRECT_DEATHS, _ageGroup, 1 );
                if (_state & Episode::EVENT_FIRST_DAY)
                    Surveys.at(_surveyPeriod)
                        .addInt( Survey::MI_FIRST_DAY_DEATHS, _ageGroup, 1 );
            }
            else if (_state & Episode::SEQUELAE)
                Surveys.at(_surveyPeriod)
                    .addInt( Survey::MI_SEQUELAE, _ageGroup, 1 );
            // Don't care about out-of-hospital recoveries
        }
    } else if (_state & Episode::SICK) {
        // Report non-malarial fever and outcomes
        Surveys.at(_surveyPeriod)
            .addInt( Survey::MI_NON_MALARIA_FEVERS, _ageGroup, 1 );

        if (_state & Episode::DIRECT_DEATH) {
            Surveys.at(_surveyPeriod)
                .addInt( Survey::MI_NMF_DEATHS, _ageGroup, 1 );
        }
    }

}

void Episode::operator& (istream& stream) {
    _time & stream;
    if (_time != TimeStep::never) {
        _surveyPeriod & stream;
        _ageGroup & stream;
        int s;
        s & stream;
        _state = Episode::State(s);
    }
}
void Episode::operator& (ostream& stream) {
    _time & stream;
    if (_time != TimeStep::never) {
        _surveyPeriod & stream;
        _ageGroup & stream;
        _state & stream;
    }
}

}
}
