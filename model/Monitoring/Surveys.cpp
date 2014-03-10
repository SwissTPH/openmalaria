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

// Must not be included _after_ boost/math/special_functions/fpclassify.hpp
#include <boost/math/nonfinite_num_facets.hpp>

#include "Monitoring/Surveys.h"
#include "Clinical/ClinicalModel.h"
#include "util/BoincWrapper.h"
#include "util/errors.h"
#include "util/CommandLine.h"
#include "interventions/InterventionManager.hpp"
#include "Simulator.h"
#include "schema/monitoring.h"

#include <gzstream/gzstream.h>
#include <fstream>
#include <stdexcept>
#include <algorithm>

namespace OM { namespace Monitoring {
    
SurveysType Surveys;

void SurveysType::init( const scnXml::Monitoring& monitoring ){
  _surveyPeriod = 0;
  _cohortOnly = false;
  
  const scnXml::Surveys::SurveyTimeSequence& survs = monitoring.getSurveys().getSurveyTime();

  _surveysTimeIntervals.reserve (survs.size() + 1);
  TimeStep last( TimeStep::never );
  for (size_t i = 0; i < survs.size(); ++i) {
      TimeStep cur(survs[i]);
    _surveysTimeIntervals.push_back( cur );
    if( last >= cur ){
        last = TimeStep::future;
    }else{
        last = cur;
    }
  }
  if( last == TimeStep::future ){
      cerr << "Warning: survey times are not listed in increasing order; will be reordered" << endl;
      sort( _surveysTimeIntervals.begin(), _surveysTimeIntervals.end() );
  }
  _surveysTimeIntervals.push_back( TimeStep::never );
  currentTimestep = _surveysTimeIntervals[0];

  Survey::init( monitoring );

  surveys.resize (_surveysTimeIntervals.size());
  if( !Simulator::isCheckpoint() ){
    for (size_t i = 0; i < surveys.size(); ++i)
        surveys[i].allocate();
  }
  current = &surveys[0];
  
  if( monitoring.getCohortOnly().present() ){
      _cohortOnly = monitoring.getCohortOnly().get();
  } else {
      // Trap potential bug in scenario design
      if( interventions::InterventionManager::cohortEnabled() ){
          throw util::xml_scenario_error( "please specify cohortOnly=\"true/false\" in monitoring element" );
      }
  }
}

void SurveysType::incrementSurveyPeriod()
{
  currentTimestep = _surveysTimeIntervals[_surveyPeriod];
  ++_surveyPeriod;
  if (_surveyPeriod >= surveys.size())
    // In this case, currentTimestep gets set to -1 so no further surveys get taken
    _surveyPeriod = 0;
  current = &surveys[_surveyPeriod];
}

void SurveysType::writeSummaryArrays ()
{
#ifdef WITHOUT_BOINC
  ofstream outputFile;		// without boinc, use plain text (for easy reading)
#else
  ogzstream outputFile;		// with, use gzip
#endif
    
    // This locale ensures uniform formatting of nans and infs on all platforms.
    locale old_locale;
    locale nfn_put_locale(old_locale, new boost::math::nonfinite_num_put<char>);
    outputFile.imbue( nfn_put_locale );

    string output_filename = util::BoincWrapper::resolveFile (util::CommandLine::getOutputName());

  outputFile.open (output_filename.c_str(), ios::out | ios::binary);

  outputFile.width (0);
  // For additional control:
  //   outputFile.precision (6);
  //   outputFile << scientific;

  for (size_t i = 1; i < surveys.size(); ++i)
    surveys[i].writeSummaryArrays (outputFile, i);

  //Infant mortality rate is a single number, therefore treated separately
  // Note: Storing a single value instead of one per reporting period is inconsistent with other
  // reporting, but I believe required for parameterisation.
  if (Survey::active[SM::allCauseIMR]) {
    outputFile << 1 << "\t" << 1 << "\t" << SM::allCauseIMR;
    outputFile << "\t" << Clinical::ClinicalModel::infantAllCauseMort() << lineEnd;
  }

  outputFile.close();
}

void SurveysType::checkpoint (istream& stream) {
    currentTimestep & stream;
    _surveysTimeIntervals & stream;
    _surveyPeriod & stream;
    // read those surveys checkpointed, call allocate on the rest:
    for( size_t i = 1; i <= _surveyPeriod; ++i )
        surveys[i] & stream;
    surveys[0].allocate();
    for( size_t i = _surveyPeriod + 1; i < surveys.size(); ++i )
        surveys[i].allocate();
    current = &surveys[_surveyPeriod];
}
void SurveysType::checkpoint (ostream& stream) {
    currentTimestep & stream;
    _surveysTimeIntervals & stream;
    _surveyPeriod & stream;
    // checkpoint only those surveys used; exclude 0 since that's a "write only DB"
    for( size_t i = 1; i <= _surveyPeriod; ++i )
        surveys[i] & stream;
}

} }
