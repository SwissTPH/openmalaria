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
#include "global.h"
#include <vector>
#include <fstream>

using namespace std;

// Forward Declaration
class Event;

//! It stores survey data and performes file I/O 
/*!
  It manages data from parasitological surey.
*/
class Summary {

 public:
  
  Summary();
  virtual ~Summary();

  void printHosts(ostream& out);

  //! Initialisation Routine. This will have to be moved to the class constructor
  void initialiseSummaries ();

  //! Clearing Routine. This will have to be moved to the class destructor
  void clearSummaryParameters();

  //! Called for all episodes which are reported. Adds an episode 
  //! to the summary data structure.
  /*  
      \param diagnosis the final diagnosis
      \param outcome the clinical outcome
      \param ageGrp age group of the patient
      \param surveyPer survey period. 
      Events are assigned to a period between 2 surveys, period 1 is
      between start of mainsimulation and the first survey. Events after the last 
    are assigned to survey period (number of surveys)+1.
  */
  void report(Event& event);
  
  //! Report a first or second line, or inpatient treatment
  void reportTreatment(int ageGroup, int regimen);

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

 private:

  //! X-dimension of summary arrays
  /*! Is numberOfSurveys+1 for all summaries */
  int _summaryDimensionX;    
  //! Y-dimension of summary arrays 
  /*! Is noOfAgeGroups for all age-specific summaries */
  int _summaryDimensionY; 
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
  //! Time intervals for all surveys specified in the XML
  vector<int> _surveysTimeIntervals; 
  int _summaryOption; //!< Binary encoded list of outputs of interest
  //! Assimilator mode
  /*! 0 for fitting give the complete output for fitting mode 1 for predictions
      give only one column to keep file sizes small.only the value column is
      recorded */
  int _assimilatorMode; 
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
    // Number of hosts transmitting to mosquitoes
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
    // Annual average kappa
    annAvgK= 26,
    // Number of episodes (non-malaria fever)
    nNMFever= 27
  };
  
  inline int isOptionIncluded (int allOptions, int option) {
    return allOptions & (1 << option);
  };
};
  
  /** Line end character. Use Unix line endings to save a little size. */
  const char lineEnd = '\n';

  template <class T>
    void writeArray(ostream& file, int measure, 
                    int assimilationMode, vector< vector<T> >& array);

  template <class T>
    void writeArray(ostream& file, int measure, 
                    int assimilationMode, vector<T>& array);



#endif
