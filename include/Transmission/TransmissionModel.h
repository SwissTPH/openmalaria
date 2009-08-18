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
  static TransmissionModel* createTransmissionModel (const std::list<Human>& population, int populationSize);
  
  //! Reads all entomological parameters from the input datafile. 
  TransmissionModel();
  //!Deallocate memory for TransmissionModel parameters and clean up
  virtual ~TransmissionModel();
  
  //FIXME: do implementing models need checkpointing?
  void write(ostream& out) const;
  void read(istream& in);
  //@}
  
  /// Set a couple of summary items
  void summarize (Summary&);
  
  
  /** Initialise the main simulation.
   *
   * Although we should have (population.size() == populationSize), it appears
   * that it's better not to use population.size(). */
  virtual void initMainSimulation (const std::list<Human>& population, int populationSize)=0; 
  
  /// Needs to be called each step of the simulation
  virtual void advancePeriod (const std::list<Human>& population, int simulationTime) {}
  
  /** Set kappa for current interval in year from infectiousness of humans.
   *
   * Also updates _annualAverageKappa.
   * 
   * NOTE: could be combined with advancePeriod(), but is currently called at
   * a different time. */
  void updateKappa (double sumWeight, double sumWt_kappa);
  
  /** Little function to copy kappa to initialKappa. */
  virtual void copyToInitialKappa () {}
  
  /** Returns the EIR, per host and per time step.
   *
   * During the pre-intervention phase, the EIR is forced, using values from
   * the XML file (possibly generated from fourier coefficients). During the
   * main simulation phase, it may be calculated or obtained from data in the
   * XML file. */
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
  
  //TODO: checkpoint
  /// Sum of all EIR returned in this timestep
  double timeStepTotalEir;
  /// Divisor of timeStepTotalEir to get average.
  int timeStepTotalEirEntries;
  
  /** Average EIR exhibited over the last year per day. */
  vector<double> eirPerDayOfYear;
};

#endif
