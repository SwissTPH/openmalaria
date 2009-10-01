/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
 OpenMalaria is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or (at
 your option) any later version.
 
 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifndef Hmod_TransmissionModel 
#define Hmod_TransmissionModel 

#include "global.h"
#include <string.h>
#include "human.h"
#include "scenario.hxx"

// Define these to print out various arrays:
//#define TransmissionModel_PrintSmoothArray

class Summary;
class PerHostTransmission;

//! Abstract base class, defines behaviour of transmission models
class TransmissionModel {
public:
  ///@brief Creation, destruction and checkpointing
  //@{
  /// Creates a derived class
  static TransmissionModel* createTransmissionModel ();
  
  //! Reads all entomological parameters from the input datafile. 
  TransmissionModel();
  //!Deallocate memory for TransmissionModel parameters and clean up
  virtual ~TransmissionModel();
  
  /** Extra initialisation, requiring information from the human population
   * structure. Only for Vector model. */
  virtual void setupNv0 (const std::list<Human>& population, int populationSize) {}
  
  void write(ostream& out) const;
  void read(istream& in);
  //@}
  
  /// Set a couple of summary items
  void summarize (Summary&);
  
  
  /** How many intervals are needed for vector initialisation before
   * updateOneLifespan() ? */
  virtual int vectorInitDuration () {
    return 0;	// Non-vector model doesn't need it
  }
  /** Called at end of vector initialisation stage. */
  virtual void endVectorInitPeriod () {}
  
  /** Initialise the main simulation.
   *
   * Although we should have (population.size() == populationSize), it appears
   * that it's better not to use population.size(). */
  virtual void initMainSimulation ()=0; 
  
  /// Needs to be called each step of the simulation
  virtual void advancePeriod (const std::list<Human>& population, int simulationTime) {}
  
  /** Set kappa for current interval in year from infectiousness of humans.
   *
   * Also updates _annualAverageKappa, and reporting of entomological
   * innoculations.
   * 
   * NOTE: could be combined with advancePeriod(), but is currently called at
   * a different time. */
  void updateKappa (double sumWeight, double sumWt_kappa);
  
  /** Little function to copy kappa to initialKappa. */
  virtual void copyToInitialKappa () {}
  
  /** Returns the EIR, per host and per time step.
   *
   * Non-vector:
   * During the pre-intervention phase, the EIR is forced, using values from
   * the XML file. During the main simulation phase, it may be calculated or
   * obtained from data in the XML file.
   *
   * Vector:
   * During vector initialisation phase, EIR is forced based on the EIR given
   * in the XML file as a Fourier Series. After endVectorInitPeriod() is called
   * the simulation switches to using dynamic EIR. advancePeriod _must_ be
   * called before this function in order to return the correct value. */
  double getEIR (int simulationTime, PerHostTransmission& host, double ageInYears);
  
  /** Set the larviciding intervention params. */
  virtual void intervLarviciding (const scnXml::Larviciding&);
  
protected:
  /** Calculates the EIR (in adults), during the main simulation phase.
   * 
   * \param simulationTime Time since start of simulation.
   * \param host The human to calculate EIR for (not used by all models). */
  virtual double calculateEIR(int simulationTime, PerHostTransmission& host, double ageInYears) = 0; 
  
  /** The type of EIR calculation. */
  int simulationMode;
  
  /** EIR per time step during the pre-intervention phase.
   *
   * Not checkpointed; doesn't need to be except when a changeEIR intervention
   * occurs. */
  vector<double> initialisationEIR; 

  /** kappa[] is the probability of infection of a mosquito at each bite.
   * It is calculated as the average infectiousness per human.
   * 
   * Checkpointed. */
  vector<double> kappa; 
  
private:
  /*! annAvgKappa is the overall proportion of mosquitoes that get infected
   * allowing for the different densities in different seasons (approximating
   * relative mosquito density with the EIR).
   *
   * Checkpointed. */
  double _annualAverageKappa;
  
  /*! Used to calculate annAvgKappa.
   *
   * Checkpointed. */
  double _sumAnnualKappa;
  
protected:
  /*! Total annual EIR.
   *
   * Checkpointed. */
  double annualEIR;
  
  /** @brief Variables for reporting of entomological innoculations to humans.
   *
   * innoculationsPer... arrays are checkpointed; timesStep... variables aren't
   * since they're calculated per step. */
  //@{
  /** Innoculations per human (all ages) per day of year.
   * 
   * Contains values from today to the previous timestep, one year ago,
   * including the initialisation phase. */
  vector<double> innoculationsPerDayOfYear;
  
  /** The total number of innoculations per age group, summed over the
   * reporting period. */
  vector<double> innoculationsPerAgeGroup;
  
  /// Sum of all EIR returned in this timestep, per age group
  vector<double> timeStepEntoInnocs;
  /// Total number of EIRs output in the timestep (roughly equal to populationSize)
  size_t timeStepNumEntoInnocs;
  //@}
};

#endif
