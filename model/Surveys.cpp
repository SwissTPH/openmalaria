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

#include "Surveys.h"
#include "inputData.h"
#include "util/BoincWrapper.h"
#include "Clinical/ClinicalModel.h"

#include "gzstream.h"
#include <fstream>
#include <stdexcept>

namespace OM {
    
SurveysType Surveys;

void SurveysType::init ()
{
  _surveyPeriod = 0;

  const scnXml::Surveys::SurveyTimeSequence& survs = InputData().getMonitoring().getSurveys().getSurveyTime();

  _surveysTimeIntervals.resize (survs.size() + 1);
  for (size_t i = 0; i < survs.size(); i++) {
    _surveysTimeIntervals[i] = survs[i];
  }
  _surveysTimeIntervals[survs.size()] = -1;
  currentTimestep = _surveysTimeIntervals[0];

  Survey::init ();

  _survey.resize (survs.size() + 1);
  for (size_t i = 0; i < _survey.size(); ++i)
    _survey[i].allocate();
  current = &_survey[0];
}

void SurveysType::incrementSurveyPeriod()
{
  currentTimestep = _surveysTimeIntervals[_surveyPeriod];
  _surveyPeriod++;
  if (_surveyPeriod >= (int) _survey.size())
    // In this case, currentTimestep gets set to -1 so no further surveys get taken
    _surveyPeriod = 0;
  current = &_survey[_surveyPeriod];
}

void SurveysType::writeSummaryArrays ()
{
#ifdef WITHOUT_BOINC
  ofstream outputFile;		// without boinc, use plain text (for easy reading)
  const char* fname = "output.txt";
#else
  ogzstream outputFile;		// with, use gzip
  const char* fname = "output.txt.gz";
#endif

  string output_filename = util::BoincWrapper::resolveFile (fname);
  ifstream test (output_filename.c_str());
  if (test.is_open())
    throw runtime_error ("File output.txt exists!");

  outputFile.open (output_filename.c_str(), ios::out | ios::binary);

  outputFile.width (0);
  // For additional control:
  //   outputFile.precision (6);
  //   outputFile << scientific;

  for (size_t i = 1; i < _survey.size(); ++i)
    _survey[i].writeSummaryArrays (outputFile, i);

  //Infant mortality rate is a single number, therefore treated separately
  // Note: Storing a single value instead of one per reporting period is inconsistent with other
  // reporting, but I believe required for parameterisation.
  if (Survey::active[imr_summary]) {
    if (!Survey::_assimilatorMode)
      outputFile << 1 << "\t" << 1 << "\t" << imr_summary;
    outputFile << "\t" << Clinical::ClinicalModel::infantAllCauseMort() << lineEnd;
  }

  outputFile.close();
}

void SurveysType::checkpoint (istream& stream) {
    currentTimestep & stream;
    _surveysTimeIntervals & stream;
    _surveyPeriod & stream;
    _survey & stream;
}
void SurveysType::checkpoint (ostream& stream) {
    currentTimestep & stream;
    _surveysTimeIntervals & stream;
    _surveyPeriod & stream;
    _survey & stream;
    current = &_survey[_surveyPeriod];
}

}