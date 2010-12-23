/* This file is part of OpenMalaria.
 *
 * Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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
#include <assert.h>

namespace OM { namespace Monitoring {
    
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
    /** Points to _survey[_surveyPeriod] (the dummy element _survey[0] before
     * start of main sim and after completion of last survey).
     * This is for data being collected for the next survey. */
    Survey* current;
    /** Timestep the current survey ends at.
     * 
     * For point-time surveys this is the time of the survey; where data is
     * collected over a period, the period is from the timestep following the
     * previous survey (or the start of the main simulation) until this time. */
    int currentTimestep;
    /** Get survey for a human -- switched depending on whether
     * the human is in the cohort and reporting mode. */
    inline Survey& getSurvey( bool inCohort ){
	if( _cohortOnly && !inCohort )
	    return _survey[0];	// output goes to dummy survey: is deleted
	else return *current;
    }
    /** As getSurvey(), except returns a number instead of a reference. */
    inline int getSurveyNumber( bool inCohort ){
        if( _cohortOnly && !inCohort )
            return 0;
        else
            return _surveyPeriod;
    }

    /** Read in some params from XML and allocate memory. */
    void init ();
    
    //! It increments the survey period
    void incrementSurveyPeriod();

    //! Write all the summary arrays requested by summaryOption to output.txt
    void writeSummaryArrays();

    ///@brief Getter functions
    //@{
    /** Return Survey number n (counting from 1).
     *
     * Survey at n=0 is the junk-survey (not reported).
     * Checks n is valid in debug mode. */
    inline Survey& at (size_t n) {
      assert (n < _survey.size());
      return _survey[n];
    }
    
    /** Return timestep of the final survey.
     *
     * We use this to control when the simulation ends.
     * This isn't quite the same as before when the simulation end was
     * explicitly specified and has a small affect on
     * infantAllCauseMortality (survey 21) output. */
    inline int getFinalTimestep () {
      return _surveysTimeIntervals[_survey.size()-2];	// final entry is a concatenated -1
    }
    //@}
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
	checkpoint (stream);
    }
    
  private:
    void checkpoint (istream& stream);
    void checkpoint (ostream& stream);
    
    //! Time intervals for all surveys specified in the XML, appended with -1
    vector<int> _surveysTimeIntervals;
    /** Index for the time dimention of the summary arrays
     * Index starts from 1 for used surveys; is 0 to write to dummy survey. */
    int _surveyPeriod;
    
    /// If true, many outputs only come from humans in the cohort.
    bool _cohortOnly;

    /// Our collection of surveys. _survey[0] is a dummy container for data
    /// we're not interested in, in order to avoid having to check current is valid.
    vector<Survey> _survey;
};
/// Data — entry-point for using Surveys. Checkpointed.
extern SurveysType Surveys;
} }
#endif
