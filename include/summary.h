/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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
#ifndef Hmod_summary
#define Hmod_summary
#include <fcntl.h>
#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Global.h"
#include <vector>
#include <fstream>

using namespace std;

// Forward Declaration
class Episode;

//! It stores survey data and performes file I/O 
/*!
  It manages data from parasitological surey.
*/
class Summary {

 public:
  
  Summary();
  virtual ~Summary();

  //! Initialisation Routine. This will have to be moved to the class constructor
  void initialiseSummaries ();

  //! Clearing Routine. This will have to be moved to the class destructor
  void clearSummaryParameters();

  /** Report a clinical episode.
   *
   * From event.getState(), an episode is reported based on severity (SICK,
   * MALARIA or COMPLICATED), and any outcomes are reported: RECOVERY (in
   * hospital, i.e. with EVENT_IN_HOSPITAL, only), SEQUELAE and DIRECT_DEATH
   * (both in and out of hospital). */
  void report(Episode& event);
  
  //! Report a first or second line, or inpatient treatment
  void reportTreatment(int ageGroup, int regimen);
  
  /// Report an indirect death (separated from report(Event) due to independant usage)
  void reportIndirectDeath (double age) {
    if (_surveyPeriod < 0) return;
    _numIndirectDeaths[_surveyPeriod][ageGroup(age)]++;
  }

  //! Report every vaccine dose given via EPI (differs from Mass for costing)
  /*!
    \param ageGroup
  */
  void reportEPIVaccination(int ageGroup);

  //! Report every vaccine dose given via campaign
  /*!
    \param ageGroup
  */
  void reportMassVaccination(int ageGroup);

  //! Report every IPT dose
  /*!
    \param ageGroup
  */
  void reportIPTDose(int ageGroup);

  //! Write all the summary arrays requested by summaryOption to output.txt
  void writeSummaryArrays();
  double infantAllCauseMort();

  //! Returns the age group/interval for a given age
  /*!
    \param age Given age
    \return age interval for a given age
  */
  int ageGroup(double age);

  //! Update(add) hosts 
  /*!
    \param age Given age
    \param value Value to add
  */
  void addToHost(double age, int value);


  //! Update(add) the infected hosts table
  /*!
    \param age Given age
    \param value Value to add
  */
  void addToInfectedHost(double age, int value);

  //! Update(add) the number of total infections
  /*!
    \param age Given age
    \param value Value to add
  */
  void addToTotalInfections(double age, int value);

  //! Update(add) the number of total patent infections
  /*!
    \param age Given age
    \param value Value to add
  */
  void addToTotalPatentInfections(double age, int value);

  //! Update(add) the table of patent hosts 
  /*!
    \param age Given age
    \param value Value to add
  */
  void addToPatentHost(double age, int value);

  //! Update(add) the table for the sum of the logarithm of the density
  /*!
    \param age Given age
    \param value Value to add
  */
  void addToSumLogDensity(double age, double value);

  //! Update(add) the table for the expected infected
  /*!
    \param age Given age
    \param value Value to add
  */
  void addToExpectedInfected(double age, double value);

  //! Update(add) the table for the pyrogenic threshold
  /*!
    \param age Given age
    \param value Value to add
  */
  void addToPyrogenicThreshold(double age, double value);

  //! Update(add) the table for the sum of log of pyrogen threshold
  /*!
    \param age Given age
    \param value Value to add
  */
  void addToSumX(double age, double value);
  
  //! It increments the survey period
  void incrementSurveyPeriod() {_surveyPeriod++;}

  //! It returns the time interval for a given survey
  /*!
    \param survey Given survey
  */
  int getSurveyTimeInterval(int survey) {return _surveysTimeIntervals[survey];};

  //! Getters & Setters
  int getNumOfAgeGroups() {return _numOfAgeGroups;};
  double getNonMalariaMortality() {return _nonMalariaMortality;};
  int getSurveyPeriod() {return _surveyPeriod;};

  void setAnnualAverageKappa(double kappa);
  void setNumTransmittingHosts(double value);
  
  void setInnoculationsPerDayOfYear (vector<double>& v) {
    _innoculationsPerDayOfYear[_surveyPeriod] = v;
  }
  void setKappaPerDayOfYear (vector<double>& v) {
    _kappaPerDayOfYear[_surveyPeriod] = v;
  }
  void setInnoculationsPerAgeGroup (vector<double>& v) {
    _innoculationsPerAgeGroup[_surveyPeriod] = v;	// copies v, not just its reference
  }
  
 private:

  //! X-dimension of summary arrays
  /*! Is numberOfSurveys+1 for all summaries */
  int _summaryDimensionX;    
  
  ///@brief All of these arrays are per survey-period (length numberOfSurveys+1, outer array), and most are per age-group (length _numOfAgeGroups).
  //@{
  vector< vector<int> > _numHosts; //!< number of hosts
  vector< vector<int> > _numInfectedHosts; //!< number of infected hosts 
  //! expected number of infected hosts //!< expected number of infected hosts 
  vector< vector<double> > _numExpectedInfected; 
  vector< vector<int> > _numPatentHosts; //!< number of patent hosts 
  vector< vector<double> > _sumX; //!< sum of log of pyrogen threshold
  //! sum of the logarithm of the density
  vector< vector<double> > _sumLogDensity; 
  vector< vector<int> > _totalInfections; //!< total infections
  //! number of hosts transmitting to mosquitoes 
  /*!  (i.e. sum of proportion of mosquitoes that
    get infected) We don't want this by age, so
    we just store kappa at surveys in a 1-dim
    array */
  //! Number of hosts transmitting to mosquitoes
  vector< double > _numTransmittingHosts;                    
  vector< vector<int> > _totalPatentInfections; //!< Total patent infections
  //!< Contribution to immunity functions
  vector< vector<double> > _contributionImmunity; 
  vector< vector<double> > _pyrogenicThreshold; //!< Sum of pyrogenic threshold
  vector< vector<int> > _numTreatments1; //!< number of treatments (1st line)
  vector< vector<int> > _numTreatments2; //!< number of treatments (2nd line)
  vector< vector<int> > _numTreatments3; //!< number of treatments (inpatient)
  //! number of episodes (uncomplicated)
  vector< vector<int> > _numUncomplicatedEpisodes; 
  vector< vector<int> > _numSevereEpisodes; //!< number of episodes (severe)
  vector< vector<int> > _numSequelae; //!< cases tith sequelae
  vector< vector<int> > _numHospitalDeaths; //!< deaths in hospital
  vector< vector<int> > _numIndirectDeaths; //!< number of deaths (indirect)
  vector< vector<int> > _numDirectDeaths; //!< number of deaths (direct)
  vector< vector<int> > _numEPIVaccines; //!< number of EPI vaccine doses
  //! number of Mass / Campaign vaccine doses
  vector< vector<int> > _numMassVaccines; 
  vector< vector<int> > _numHospitalRecoveries; //!< recoveries in hospital
  vector< vector<int> > _numHospitalSequelae; //!< sequelae in hospital
  vector< vector<int> > _numIPTDoses; //!< number of IPT Doses
  vector< double > _annualAverageKappa; //!< Annual Average Kappa
  //! Number of episodes (non-malaria fever)
  vector< vector<int> > _numNonMalariaFever; 
  /// Innoculations per human (all ages) per day of year, over the last year.
  vector< vector<double> > _innoculationsPerDayOfYear;
  /// Kappa (human infectiousness) weighted by availability per day-of-year for the last year.
  vector< vector<double> > _kappaPerDayOfYear;
  /** The total number of innoculations per age group, summed over the
   * reporting period. */
  vector< vector<double> > _innoculationsPerAgeGroup;
  //@}
  
  //! Time intervals for all surveys specified in the XML
  vector<int> _surveysTimeIntervals; 
  int _summaryOption; //!< Binary encoded list of outputs of interest
  /** Assimilator mode
   *
   * If true, skip the first 3 columns of output to reduce file size. */
  bool _assimilatorMode; 
  //! Number of agegroups given in the XML
  //! Plus 1 for individual older than the highest upperbound
  int _numOfAgeGroups; 
  double _lowerbound; //!< Lower boundary of the youngest agegroup
  //! Upper boundary of agegroups, in years.
  //! Upperbound is lowerbound for the next older agegroup
  vector<double> _upperbound; 
  //! Index for the time dimention of the summary arrays
  /*! To store events that happen between to surveys.Survey is always the
      (one-based) index of the previous survey plus one, ie. 1 before the first
      survey, numberOfSurveys+1 for the time after the last survey.
  */
  int _surveyPeriod; 
  double _nonMalariaMortality; //!< Non-malaria mortality in under 1year olds.
  
  /// Enumeration of reporting options
  enum Measure{

    nHost = 0,
    nInfect = 1,
    // expected number of infected hosts
    nExpectd= 2,
    // number of patent hosts
    nPatent= 3,
    // sum of log of pyrogen threshold
    sumX= 4,
    // sum of the logarithm of the density
    sumlogDens= 5,
    // Total infections
    totalInfs= 6,
    /// Infectiousness of human population to mosquitoes
    nTransmit= 7,
    // total patent infections
    totalPatentInf= 8,
    // contribution to immunity functions
    contrib= 9,
    // sum of pyrogenic threshold		
    pyrogenThrs= 10,
    // number of treatments (1st line)
    nTreatments1= 11,
    // number of treatments (2nd line)
    nTreatments2= 12,
    // number of treatments (inpatient)
    nTreatments3= 13,
    // number of episodes (uncomplicated)
    nUncomp= 14,
    // number of episodes (severe)
    nSevere= 15,
    // cases with sequelae      
    nSeq= 16,
    // deaths in hospital
    nHospitalDeaths= 17,
    // number of deaths (indirect)
    nIndDeaths= 18,
    // number of deaths (direct)
    nDirDeaths= 19,
    // number of EPI vaccine doses
    nEPIVaccines= 20,
    //all cause infant mortality rate
    imr_summary= 21,
    // number of Mass/Campaign vaccine doses
    nMassVaccines= 22,
    // recoveries in hospital
    nHospitalRecovs= 23,
    // sequelae in hospital
    nHospitalSeqs= 24,
    // IPT doses
    nIPTDoses= 25,
    /// Calculated once a year as sum of human infectiousness divided by initial EIR summed over a year
    annAvgK= 26,
    // Number of episodes (non-malaria fever)
    nNMFever= 27,
    // EIR per day of year, summed over all years
    innoculationsPerDayOfYear = 28,
    // Kappa per day of year, for the last year
    kappaPerDayOfYear = 29,
    /** Report the total number of innoculations per age group, summed over the
     * reporting period. */
    innoculationsPerAgeGroup = 30,
    
    // Note: can't use values greater than 31 without forcing a 64-bit type
  };
  
  inline int isOptionIncluded (int allOptions, int option) {
    return allOptions & (1 << option);
  };
};

/** Line end character. Use Unix line endings to save a little size. */
const char lineEnd = '\n';

template <class T>
void writeArray(ostream& file, int measure, bool assimilationMode, vector< vector<T> >& array);

template <class T>
void writeArray(ostream& file, int measure, bool assimilationMode, vector<T>& array);

#endif
