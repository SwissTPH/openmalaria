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
#include "Clinical/CaseManagementCommon.h"
#include "Host/Human.h"

namespace OM {
namespace Clinical {
using namespace Monitoring;


Episode::~Episode ()
{
    report ();
}

void Episode::flush() {
    report();
    time = sim::never();
}


void Episode::update (const Host::Human& human, Episode::State newState)
{
    if( time + healthSystemMemory < sim::now() ){
        report ();

        time = sim::now();
        surveyPeriod = Survey::getSurveyNumber();
        ageGroup = human.getMonitoringAgeGroup();
        cohortSet = human.cohortSet();
        state = newState;
    } else {
        state = Episode::State (state | newState);
    }
}

void Episode::report () {
    if (time < sim::zero())        // Nothing to report
        return;
    
    Survey& survey = Survey::getSurvey(surveyPeriod);

    // Reports malarial/non-malarial UC fever dependent on cause, not diagnosis.
    if (state & Episode::MALARIA) {
        // Malarial fevers: report bout
        if (state & Episode::COMPLICATED) {
            survey.addInt( Report::MI_SEVERE_EPISODES, ageGroup, cohortSet, 1) ;
        } else { // UC or UC2
            survey.addInt( Report::MI_UNCOMPLICATED_EPISODES, ageGroup, cohortSet, 1 );
        }

        // Report outcomes of malarial fevers
        if (state & Episode::EVENT_IN_HOSPITAL) {
            if (state & Episode::DIRECT_DEATH) {
                survey
                    .addInt( Report::MI_DIRECT_DEATHS, ageGroup, cohortSet, 1 )
                    .addInt( Report::MI_HOSPITAL_DEATHS, ageGroup, cohortSet, 1 );
                if (state & Episode::EVENT_FIRST_DAY)
                    survey
                        .addInt( Report::MI_FIRST_DAY_DEATHS, ageGroup, cohortSet, 1 )
                        .addInt( Report::MI_HOSPITAL_FIRST_DAY_DEATHS, ageGroup, cohortSet, 1 );
            }
            else if (state & Episode::SEQUELAE) {
                survey
                    .addInt( Report::MI_SEQUELAE, ageGroup, cohortSet, 1 )
                    .addInt( Report::MI_HOSPITAL_SEQUELAE, ageGroup, cohortSet, 1 );
            }
            else if (state & Episode::RECOVERY)
                survey
                    .addInt( Report::MI_HOSPITAL_RECOVERIES, ageGroup, cohortSet, 1 );
        } else {
            if (state & Episode::DIRECT_DEATH) {
                survey
                    .addInt( Report::MI_DIRECT_DEATHS, ageGroup, cohortSet, 1 );
                if (state & Episode::EVENT_FIRST_DAY)
                    survey
                        .addInt( Report::MI_FIRST_DAY_DEATHS, ageGroup, cohortSet, 1 );
            }
            else if (state & Episode::SEQUELAE)
                survey
                    .addInt( Report::MI_SEQUELAE, ageGroup, cohortSet, 1 );
            // Don't care about out-of-hospital recoveries
        }
    } else if (state & Episode::SICK) {
        // Report non-malarial fever and outcomes
        survey
            .addInt( Report::MI_NON_MALARIA_FEVERS, ageGroup, cohortSet, 1 );

        if (state & Episode::DIRECT_DEATH) {
            survey
                .addInt( Report::MI_NMF_DEATHS, ageGroup, cohortSet, 1 );
        }
    }

}

void Episode::operator& (istream& stream) {
    time & stream;
    if (time > sim::zero()) {
        surveyPeriod & stream;
        ageGroup & stream;
        cohortSet & stream;
        int s;
        s & stream;
        state = Episode::State(s);
    }
}
void Episode::operator& (ostream& stream) {
    time & stream;
    if (time >= sim::zero()) {
        surveyPeriod & stream;
        ageGroup & stream;
        cohortSet & stream;
        state & stream;
    }
}

}
}
