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
#ifndef Hmod_human
#define Hmod_human
#include "global.h"
#include "event.h"
#include "WithinHostModel.h"
#include "Pathogenesis/PathogenesisModel.h"
#include "drug.h"
#include "TransmissionModel/PerHost.h"

// Forward declaration
class TransmissionModel;
class CaseManagementModel;

//! Model of a human individual 
class Human {
public:
  
  /// @name Constructors
  //@{
  /** Initialise all variables of a human datatype.
   * 
   * \param ID unique identifier
   * \param dateOfBirth date of birth in time steps */
  Human(TransmissionModel& tm, int ID, int dateOfBirth, int simulationTime);

  /**  Initialise all variables of a human datatype including infectionlist and
   * and druglist.
   * \param funit IO unit */
  Human(istream& funit, TransmissionModel& tm, int simulationTime);
  //@}
  
  /** Destructor
   * 
   * NOTE: this destructor does nothing to allow shallow copying to the
   * population list. Human::destroy() does the real freeing and must be
   * called explicitly. */
  ~Human() {}
  
  /// The real destructor
  void destroy();
  
  /// @name Checkpointing functions
  //@{
  friend ostream& operator<<(ostream& out, const Human& human);
  //@}
  
  /** @name External functions
   * @brief Functions not called from within Human but calling functions in Human */
  //@{
  /** If the individual is too old, returns true. Otherwise, updates the
   * individual for the time-step. */
  bool update(int simulationTime, TransmissionModel* transmissionModel);
  
  //! Summarize the state of a human individual.
  void summarize();
  //@}
  
  /// @name updateInfection related functions
  //@{
  void updateInfection(double expectedInfectionRate, double expectedNumberOfInfections);
  //@}
  
  /// @name determineClinicalStatus related functions
  //@{
  void determineClinicalStatus();
  //@}
  
  /// @name updateInterventionStatus related functions
  //@{
  /*! Apply interventions to this human if eligible. Calculate the remaining
      efficacy of the latest vaccination if vaccinated before */
  void updateInterventionStatus();
  
  /*! 
    Update the number of doses and the date of the most recent vaccination in
    this human */
  void vaccinate();
  //@}
  
  /** @name More functions
   * @brief Functions used internally by more than one category above
   * (excluding update and summarize) */
  //@{
  //! Determines the age group of a human
  int ageGroup() const;
  
  //! Get the age in years, based on an input reference time.
  double getAgeInYears(int time) const;

  //! Get the age in years, based on current simulationTime in human.
  double getAgeInYears() const;
  //@}
  
  /** @name Other functions
   * @brief Other functions not called within Human */
  //@{
  //! Get the PEV Efficacy
  double getPEVEfficacy() {return _PEVEfficacy;} 
  
  //! Get the Availability to mosquitoes
  double getBaselineAvailabilityToMosquitoes() const {return _perHostTransmission._BaselineAvailabilityToMosquitoes;}
  
  inline double getProbTransmissionToMosquito() const {return _withinHostModel->getProbTransmissionToMosquito();}
  
  //! Return doomed value
  int getDoomed() {return _doomed;};

  //! Returns the date of birth
  int getDateOfBirth() {return _dateOfBirth;};
  
  double getTotalDensity() const {return _totalDensity;}
  
  void IPTiTreatment (double compliance);
  
  void setProbabilityOfInfection(double probability) { _perHostTransmission._pinfected=probability;};
  //@}
  
  /** @name OldWithinHostModel functions
   * @brief Functions only used by oldWithinHostModel.cpp */
  //@{
  double getTimeStepMaxDensity() const {return _timeStepMaxDensity;}
  void setTimeStepMaxDensity(double timeStepMaxDensity){_timeStepMaxDensity = timeStepMaxDensity;}
  
  void setTotalDensity(double totalDensity) {_totalDensity=totalDensity;}
  
  double getBSVEfficacy() {return _BSVEfficacy;}
  
  /*!  Determines the probability that the individual transmits to a feeding
    mosquito */
  double infectiousness();

  //@}
  
  /// For direct interactions with within host model
  ///TODO: possibly functionality of functions using them should be moved to Human?
  //@{
  int getCumulativeInfections() {return _withinHostModel->getCumulativeInfections();}
  void clearInfections();
  //@}
  
  ///@name static public
  //@{
  static void initHumanParameters ();
  
  static void clear();
  
/*
  The detection limit (in parasites/ul) is currently the same for PCR and for microscopy
  TODO: in fact the detection limit in Garki should be the same as the PCR detection limit
  The density bias allows the detection limit for microscopy to be higher for other sites
*/
  static double detectionlimit;
  //@}
  
  PerHostTransmission _perHostTransmission;
  
private:
  ///@name Private variables
  //@{
  // Time from start of the simulation
  int _simulationTime;
  
  WithinHostModel* _withinHostModel;

  CaseManagementModel* _caseManagement;
  
  PathogenesisModel* _pathogenesisModel;
  
  //!Total asexual blood stage density
  double _ylag[4];
  //!Date of birth, time step since start of warmup
  int _dateOfBirth;
  //!Indicates that individual will die from indirect mortality
  int _doomed;
  //!unique identifier
  int _ID;
  //!number of vaccine doses this individual has received
  int _lastVaccineDose;
  //!Total asexual blood stage density
  double _totalDensity;		// possibly move to WithinHostModel
  //!Remaining efficacy of Blood-stage vaccines
  double _BSVEfficacy;
  //!Remaining efficacy of Pre-erythrocytic vaccines
  double _PEVEfficacy;
  //!Remaining efficacy of Transmission-blocking vaccines
  double _TBVEfficacy;
  //!Maximum parasite density during the previous 5-day interval
  double _timeStepMaxDensity;	// WithinHostModel, used by PathogenesisModel
  //@}
  
  //! Determines eligibility and gives IPTi SP or placebo doses 
  void setLastSPDose();

  void updateInfectionStatus();

  //! docu 
  void computeFinalConc(Drug drg, int hl);

  //! In practice calculates drug concs (PK) for each drug
  void startWHCycle();

  double convertDose(Drug drg, Dose ds);


  //! docu
  double computeExponentialDecay(double c, int hl, int t);

  void clearInfection(Infection *iCurrent);

  //! Reads drugs from checkpoint
  void readDrugs(fstream& funit, int multiplicity);

  //!  Write all drugs of this human to the standard checkpoint file.
  /*
    \param funit io unit number
  */
  void writeDrugs(fstream& funit);
};

#endif
