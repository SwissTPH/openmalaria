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

#ifndef Hmod_Surveys
#define Hmod_Surveys

#include "Monitoring/Survey.h"
#include "Global.h"
#include "Host/Human.h"
#include <assert.h>

namespace OM {
namespace interventions {
    class InterventionManager;
}
namespace Monitoring {
    
/** Class to collect surveys and write them out.
 *
 * Surveys are written to the file output.txt. There is a one-to-one mapping of
 * lines to data entries (except the file ends with a new line). Data columns are
 * separated by tabs on each.
 * 
 * The first column lists the survey number as an integer, counting from one,
 * the second column a "group" parameter as a string (precise meaning depends
 * on the measure), the third column the measure as an integer ID (the values
 * in the SurveyCodes enum), and the forth a value (integer or floating-point,
 * but when exported to the database always considered a double). */
class SurveysType
{
public:
    SurveysType() :
            currentTimestep(TimeStep::never),
            m_cohortOnly(false)
            {}
    
    /** Timestep the current survey ends at.
     * 
     * For point-time surveys this is the time of the survey; where data is
     * collected over a period, the period is from the timestep following the
     * previous survey (or the start of the main simulation) until this time. */
    TimeStep currentTimestep;
    inline bool getCohortOnly()const{ return m_cohortOnly; }

    /** Read in some params from XML and allocate memory. */
    void init (const scnXml::Monitoring& monitoring);
    
    //! It increments the survey period
    void incrementSurveyPeriod();

    //! Write all the summary arrays requested by summaryOption to output.txt
    void writeSummaryArrays();

    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        checkpoint (stream);
    }
    
private:
    void checkpoint (istream& stream);
    void checkpoint (ostream& stream);
    
    //! Time intervals for all surveys specified in the XML, appended with -1
    vector<TimeStep> _surveysTimeIntervals;
    
    /// If true, many outputs only come from humans in the cohort.
    bool m_cohortOnly;

    /// Our collection of surveys. m_surveys[0] is a dummy container for data
    /// we're not interested in, in order to avoid having to check current is valid.
    vector<Survey> m_surveys;
    
    friend class Survey;
};
/// Data — entry-point for using Surveys. Checkpointed.
extern SurveysType Surveys;
} }
#endif
