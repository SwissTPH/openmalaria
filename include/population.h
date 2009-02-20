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
#ifndef Hmod_population
#define Hmod_population
#include "global.h"
#include <list>
#include <fstream>
#include <vector>

/*
  This is the maximum age of an individual that the simulation program can handle.
  Max age for a scenario is given in the  the xml file.
*/
const int maxLifetimeDays= 32855;
const int ngroups= 20;

//a0 and a1 are the lower and upper age limits for the field data demography agegroups
extern double a0[ngroups];
extern double a1[ngroups];
extern double perc[ngroups];
//demography variables used in estimating the smooth curve
extern double M1[ngroups];
extern double M2[ngroups];
extern double M[ngroups];
extern double pred[ngroups];
//parameters defining smooth curve of target age-distribution
extern double mu0;
extern double mu1;
extern double alpha0;
extern double alpha1;
extern double rho;
/*
  target cumulative percentage of population by age
  TODO5D
*/
extern double *cumpc;
extern int cumpcX;

extern int IDCounter; 

// Forward declarations
class Human;
class CaseManagement;
class TransmissionModel;

//! The simulated human population
class Population{

 public:

  Population();
  Population(TransmissionModel* transmissionModel, int _populationSize);

  ~Population();

   /*! Estimates demography parameters to define a smooth curve for the target
      population age-distribution (age in years) */
  void estimateRemovalRates();

  /*! Takes the best-fitting demography parameters estimated by
      estimateRemovalRates and sets up the initial population according to
      these */
  void setupPyramid(bool isCheckpoint);

  //! Initialise arrays for recording infant mortality
  void initialiseInfantArrays();
 
  //! Write the human list and some parameters to a file. This is called during checkpoints
  void writeLists(fstream& funit);

  //! Read all lists from a checkpoint file after resuming the application
  void readLists(fstream& funit);
  
  //! Updates all individuals in the list for one time-step
  /*!  Also updates the population-level measures such as infectiousness, and
       the age-distribution by c outmigrating or creating new births if
       necessary */
  void update1();
  
  //! Checks for time-based interventions and implements them
  /*!   
       \param time Current time (in tsteps) 
  */
  void implementIntervention(int time);

  //!  Mass IPTi Treatment Intervention
  /*!   
       \param time Current time (in tsteps) 
  */
  void massIPTiTreatment(int time);

  //!   Clear human collection.
  void clear();

  //! Mass Vaccination Intervention
  /*!
     \param time Current time (in tsteps)
  */
  void vaccinatePopulation(int time);
  
  //! Initialise human list
  void initialiseHumanList();
  
  //! Updates infant mortality array when there is a death
  void updateInfantArrays(int agetstep, int doomed);
 
  //! Creates initializes and add to the population list a new uninfected human
  /*! 
     \param dob date of birth (usually current time)
  */
  void newHuman(int dob);

  //! Makes a survey
  void newSurvey();
 
  //! Mass Drug Treatment Intervention
  /*!  
    \param time Current time (in tsteps)
  */
  void massTreatment(int time);

  //! remove human from the list
  short outmigrate(Human *current, int Nsize, int *survivsSoFar, int *nCounter, int *pCounter, int *k, double yage);
 
  //! estimates parameters describing smooth demography target population
  double setDemoParameters (double param1, double param2);

 private:

  //! Size of the human population   
  int _populationSize;

  //! TransmissionModel model
  TransmissionModel* _transmissionModel;

  //! Case Management System model
  CaseManagement* _caseManagement;
  
  //The simulated human population
  std::list<Human*> _population;

  //! max lifespan in intervals
  int _maxTimestepsPerLife;

  /*!
    this is needed to prevent checkpoint cheats. Ideally a unique identifier per workunit, but a
    random integer number should do the job
  */
  int _workUnitIdentifier;
  
  //! array for stored prevalences 20-25 years for 5 months (for neonatal deaths)
  std::vector<double> matArray;

  /*!
    annAvgKappa is the overall proportion of mosquitoes that get infected allowing for the different densities
    in different seasons (approximating relative mosquito density with the EIR)
  */
  double _annualAverageKappa;
 
  //! Used to calculate annAvgKappa.
  double _sumAnnualKappa;

  void init();

};

#endif
