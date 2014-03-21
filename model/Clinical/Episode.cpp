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
#include "Host/Human.h"

namespace OM {
namespace Clinical {
using namespace Monitoring;

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


void Episode::update (const Host::Human& human, Episode::State newState)
{
    if (TimeStep::simulation > (_time + healthSystemMemory)) {
        report ();

        _time = TimeStep::simulation;
        _surveyPeriod = Survey::getSurveyNumber();
        _ageGroup = human.getMonitoringAgeGroup();
        _inCohort = human.isInAnyCohort();
        _state = newState;
    } else {
        _state = Episode::State (_state | newState);
    }
}

void Episode::report () {
    if (_time == TimeStep::never)        // Nothing to report
        return;
    
    Survey& survey = Survey::getSurvey(_surveyPeriod);

    // Reports malarial/non-malarial UC fever dependent on cause, not diagnosis.
    if (_state & Episode::MALARIA) {
        // Malarial fevers: report bout
        if (_state & Episode::COMPLICATED) {
            survey.addInt( Report::MI_SEVERE_EPISODES, _ageGroup, _inCohort, 1) ;
        } else { // UC or UC2
            survey.addInt( Report::MI_UNCOMPLICATED_EPISODES, _ageGroup, _inCohort, 1 );
        }

        // Report outcomes of malarial fevers
        if (_state & Episode::EVENT_IN_HOSPITAL) {
            if (_state & Episode::DIRECT_DEATH) {
                survey
                    .addInt( Report::MI_DIRECT_DEATHS, _ageGroup, _inCohort, 1 )
                    .addInt( Report::MI_HOSPITAL_DEATHS, _ageGroup, _inCohort, 1 );
                if (_state & Episode::EVENT_FIRST_DAY)
                    survey
                        .addInt( Report::MI_FIRST_DAY_DEATHS, _ageGroup, _inCohort, 1 )
                        .addInt( Report::MI_HOSPITAL_FIRST_DAY_DEATHS, _ageGroup, _inCohort, 1 );
            }
            else if (_state & Episode::SEQUELAE) {
                survey
                    .addInt( Report::MI_SEQUELAE, _ageGroup, _inCohort, 1 )
                    .addInt( Report::MI_HOSPITAL_SEQUELAE, _ageGroup, _inCohort, 1 );
            }
            else if (_state & Episode::RECOVERY)
                survey
                    .addInt( Report::MI_HOSPITAL_RECOVERIES, _ageGroup, _inCohort, 1 );
        } else {
            if (_state & Episode::DIRECT_DEATH) {
                survey
                    .addInt( Report::MI_DIRECT_DEATHS, _ageGroup, _inCohort, 1 );
                if (_state & Episode::EVENT_FIRST_DAY)
                    survey
                        .addInt( Report::MI_FIRST_DAY_DEATHS, _ageGroup, _inCohort, 1 );
            }
            else if (_state & Episode::SEQUELAE)
                survey
                    .addInt( Report::MI_SEQUELAE, _ageGroup, _inCohort, 1 );
            // Don't care about out-of-hospital recoveries
        }
    } else if (_state & Episode::SICK) {
        // Report non-malarial fever and outcomes
        survey
            .addInt( Report::MI_NON_MALARIA_FEVERS, _ageGroup, _inCohort, 1 );

        if (_state & Episode::DIRECT_DEATH) {
            survey
                .addInt( Report::MI_NMF_DEATHS, _ageGroup, _inCohort, 1 );
        }
    }

}

void Episode::operator& (istream& stream) {
    _time & stream;
    if (_time != TimeStep::never) {
        _surveyPeriod & stream;
        _ageGroup & stream;
        _inCohort & stream;
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
        _inCohort & stream;
        _state & stream;
    }
}

}
}
