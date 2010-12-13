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

#include "Global.h"
#include "Host/Human.h"
#include "util/errors.h"
#include "scenario.hxx"

#include <fstream>
#include <string.h>

// Define these to print out various arrays:
//#define TransmissionModel_PrintSmoothArray

namespace OM {
    class Summary;
namespace Transmission {
    class PerHostTransmission;


/** There are 3 simulation modes. */
enum SimulationMode {
  /** Equilibrium mode
   * 
   * This is used for the warm-up period and if we want to separate direct
   * effect of an intervention from indirect effects via transmission
   * intensity. The seasonal pattern and intensity of the EIR do not change
   * over years.
   * 
   * For the vector model, this runs most calculations dynamically but still
   * forces the EIR. */
  equilibriumMode = 2,
  
  /** Transient EIR known
   * 
   * This is used to simulate an intervention that changes EIR, and where we
   * have measurements of the EIR over time during the intervention period. */
  transientEIRknown = 3,
  
  /** EIR changes
   * 
   * The simulation is driven by the EIR which changes dynamically during the
   * intervention phase as a function of the characteristics of the
   * interventions.
   * 
   * Dependending on whether the Vector or NonVector model is in use, this EIR
   * may be calculated from a mosquito emergence rate or be an input EIR
   * scaled by the relative infectiousness of the humans. */
  dynamicEIR = 4,
};


//! Abstract base class, defines behaviour of transmission models
class TransmissionModel {
public:
  ///@brief Creation, destruction and checkpointing
  //@{
  /// Creates a derived class
  static TransmissionModel* createTransmissionModel (int populationSize);
  
protected:
  //! Reads all entomological parameters from the input datafile. 
  TransmissionModel();
public:
  //!Deallocate memory for TransmissionModel parameters and clean up
  virtual ~TransmissionModel();
  
  /** Extra initialisation when not loading from a checkpoint, requiring
   * information from the human population structure. */
  virtual void setupNv0 (const std::list<Host::Human>& population, int populationSize) {}
  
  /// Checkpointing
  template<class S>
  void operator& (S& stream) {
      checkpoint (stream);
  }
  //@}
  
  /** Set some summary items.
   *
   * Overriding functions should call this base version too. */
  virtual void summarize (Monitoring::Survey& survey);
  
  /** Scale the EIR used by the model.
   *
   * EIR is scaled in memory (so will affect this simulation).
   * XML data is not touched. */
  virtual void scaleEIR (double factor) =0;
  
  /** Scale the EIR descriptions in the XML element.
   * This updates the XML, and not the EIR descriptions used for simulations.
   * In order for changes to be written back to the XML file,
   * InputData.documentChanged needs to be set.
   * 
   * @param ed	Access to XML element to update.
   * @param factor	Multiplicative factor by which to scale EIR. */
  virtual void scaleXML_EIR (scnXml::EntoData& ed, double factor) const =0;
  
  /** How many intervals are needed for vector initialisation after one-
   * lifespan initialisation?
   *
   * Units: timesteps; value should be a whole number of years. */
  virtual int transmissionInitDuration () =0;
  /** Called after end of transmissionInitDuration() to control initialisation
   * iterations.
   * 
   * Should determine whether another init iteration is needed, make necessary
   * adjustments, and return number of timesteps to run this initialisation for
   * (0 if no further iteration is needed).
   *
   * Units: timesteps; value should be a whole number of years. */
  virtual int transmissionInitIterate () {
    return 0;
  }
  
  /** Initialise the main simulation.
   *
   * Although we should have (population.size() == populationSize), it appears
   * that it's better not to use population.size(). */
  virtual void initMainSimulation ()=0; 
  
  /** Needs to be called each step of the simulation before nearly anything
   * else.
   * 
   * Calculates the ageCorrectionFactor which is used in much of the
   * transmission code, including getEIR (used in Human update). */
  void updateAgeCorrectionFactor (std::list<Host::Human>& population, int populationSize);
  
  /** Needs to be called each step of the simulation
   *
   * Runs internal calculations of Vector model. */
  virtual void vectorUpdate (const std::list<Host::Human>& population, int simulationTime) {};
  
  /** Needs to be called each step of the simulation
   *
   * Summarises infectiousness of humans and of mosquitoes. */
  void updateKappa (const std::list<Host::Human>& population, int simulationTime);
  
  virtual void changeEIRIntervention (const scnXml::NonVector&) {
      throw util::xml_scenario_error("changeEIR intervention can only be used with NonVectorTransmission model!");
  }
  
  /** Returns the EIR (inoculation rate per host per time step).
   *
   * Non-vector:
   * During the pre-intervention phase, the EIR is forced, using values from
   * the XML file. During the main simulation phase, it may be calculated or
   * obtained from data in the XML file.
   *
   * Vector:
   * During vector initialisation phase, EIR is forced based on the EIR given
   * in the XML file as a Fourier Series. After endVectorInitPeriod() is called
   * the simulation switches to using dynamic EIR. advanceStep _must_ be
   * called before this function in order to return the correct value. */
  double getEIR (int simulationTime, PerHostTransmission& host, double ageYears, Monitoring::AgeGroup ageGroup);
  
  /** Set the larviciding intervention params. */
  virtual void intervLarviciding (const scnXml::Larviciding&);
  
  /** Remove all current infections to mosquitoes, such that without re-
   * infection, humans will then be exposed to zero EIR. */
  virtual void uninfectVectors() =0;
  
protected:
  /** Calculates the EIR individuals are exposed to.
   * 
   * @param simulationTime Time since start of simulation.
   * @param host Transmission data for the human to calculate EIR for.
   * @param ageGroupData Age group of this host for availablility data.
   *
   * @returns  The age- and heterogeneity-specific EIR an individual is exposed
   * to, in units of inoculations per day. */
  virtual double calculateEIR(int simulationTime, PerHostTransmission& host, double ageYears) = 0; 
  
  virtual void checkpoint (istream& stream);
  virtual void checkpoint (ostream& stream);
  
private:
    void ctsCbInputEIR (ostream& stream);
    void ctsCbSimulatedEIR (ostream& stream);
    void ctsCbKappa (ostream& stream);
    void ctsCbHumanAvail (ostream& stream);
    void ctsCbNumTransmittingHumans (ostream& stream);
    
public:
  /** Correction factor for PerHostTransmission::getRelativeAvailability.
   *
   * This parameter is recalculated every time-step; dependant on the population,
   * thus it is better to store it in a per-population object than as a static
   * parameter. */
  double ageCorrectionFactor;
  
protected:
  /** The type of EIR calculation. Checkpointed. */
  int simulationMode;
  
  /** Entomological inoculation rate for adults during the
   * pre-intervention phase.
   * 
   * Length: Global::intervalsPerYear
   *
   * Units: infectious bites per adult per timestep
   *
   * Not checkpointed; doesn't need to be except when a changeEIR intervention
   * occurs. */
  vector<double> initialisationEIR; 

  /** kappa[] is the probability of infection of a mosquito at each bite.
   * It is calculated as the average infectiousness per human.
   * 
   * Checkpointed. */
  //TODO: kappa should be moved to NonVectorTransmission class. Reporting should
  // be model specific (isn't currently correct for VectorTransmission).
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
  
  double BSSInitialisationEIR;
  double BSSInnoculationsPerDayOfYear;
  int BSSTimesteps;
  
  /// For "num transmitting humans" cts output.
  int numTransmittingHumans;

protected:
  /** Total annual infectious bites per adult.
   *
   * Checkpointed. */
  double annualEIR;
  
  /** @brief Variables for reporting of entomological inoculations to humans.
   *
   * innoculationsPer... arrays are checkpointed; timesStep... variables aren't
   * since they're calculated per step. */
  //@{
  /** Inoculations per human (all ages) per day of year.
   * 
   * Contains values from today to the previous timestep, one year ago,
   * including the initialisation phase. */
  vector<double> innoculationsPerDayOfYear;
  
  /** The total number of inoculations per age group, summed over the
   * reporting period. */
  vector<double> innoculationsPerAgeGroup;
  
  /** Sum of all EIR returned in this timestep, per age group
   * Doesn't need to be checkpointed. */
  vector<double> timeStepEntoInnocs;
  /** Total number of EIRs output in the timestep (roughly equal to populationSize)
   * Doesn't need to be checkpointed. */
  size_t timeStepNumEntoInnocs;
  //@}
  
  /** @brief Variables for shared graphics kappa-by-age graph
   * Don't need checkpointing; only kept here to save reallocating each step. */
  //@{
  size_t noOfAgeGroupsSharedMem;
  vector<double> kappaByAge;
  vector<int> nByAge;
  //@}
};

} }
#endif
