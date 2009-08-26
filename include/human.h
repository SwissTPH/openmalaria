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
#include "Transmission/PerHost.h"
#include "InfectionIncidenceModel.h"
#include "WithinHost/WithinHostModel.h"

// Forward declaration
class TransmissionModel;
class ClinicalModel;

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
  Human(istream& funit, TransmissionModel& tm);
  //@}
  
  /** Destructor
   * 
   * Note: this destructor does nothing in order to allow shallow copying of a
   * Human into the population list. Human::destroy() does the real freeing and
   * must be called explicitly. */
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
  void updateInfection(TransmissionModel*);
  //@}
  
  /// @name updateInterventionStatus related functions
  //@{
  /*! Apply interventions to this human if eligible. Calculate the remaining
      efficacy of the latest vaccination if vaccinated before */
  void updateInterventionStatus();
  
  /** A wrapper around vaccinate to also report.
   * Can't move reporting to vaccinate as it would change existing reporting. */
  void massVaccinate ();
  //@}
  
  /** @name More functions
   * @brief Functions used internally by more than one category above
   * (excluding update and summarize) */
  //@{
  //! Determines the age group of a human
  int ageGroup() const;
  
  //! Get the age in years, based on current simulationTime.
  double getAgeInYears() const;
  //@}
  
  /** @name Other functions
   * @brief Other functions not called within Human */
  //@{
  //! Returns the date of birth
  int getDateOfBirth() {return _dateOfBirth;};
  
  void IPTiTreatment ();
  //@}
  
  /** @name OldWithinHostModel functions
   * @brief Functions only used by oldWithinHostModel.cpp */
  //@{
  double getBSVEfficacy() {return _BSVEfficacy;}
  
  /*!  Determines the probability that the individual transmits to a feeding
    mosquito */
  double infectiousness();

  //@}
  
  /// For direct interactions with within host model
  void clearInfections();
  
  /// Give human a new ITN
  void setupITN ();
  /// Give human a new IRS
  void setupIRS ();
  /// Give human a new VA intervention
  void setupVA ();
  
  ///@name static public
  //@{
  static void initHumanParameters ();
  
  static void clear();
  //@}
  
  
  /** @brief Models
   *
   * These contain various sub-models used by Humans. */
  //@{
  /// Contains per-species vector data (VectorTransmission only).
  PerHostTransmission perHostTransmission;
  
  /// The InfectionIncidenceModel translates per-host EIR into new infections
  InfectionIncidenceModel *infIncidence;
  
  /// The WithinHostModel models parasite density and immunity
  WithinHostModel *withinHostModel;
  
private:
  /*! Update the number of doses and the date of the most recent vaccination in
   * this human */
  void vaccinate();
  
  /** The ClinicalModel encapsulates pathogenesis (sickness status),
   * case management (medicating drugs)
   * and clinical outcomes (morbidity, reporting). */
  ClinicalModel *clinicalModel;
  //@}
  
  ///@brief Private variables
  //@{
  //!Total asexual blood stage density
  double _ylag[4];
  //!Date of birth, time step since start of warmup
  int _dateOfBirth;
  //!unique identifier
  int _ID;
  /** Number of vaccine doses this individual has received.
   *
   * If an individual misses one EPI (continuous) vaccine dose, it's
   * intentional that they also miss following EPI doses (unless a timed mass
   * vaccination reintroduces them to the EPI schedule). */
  int _lastVaccineDose;
  //!Remaining efficacy of Blood-stage vaccines
  double _BSVEfficacy;
  //!Remaining efficacy of Pre-erythrocytic vaccines
  double _PEVEfficacy;
  //!Remaining efficacy of Transmission-blocking vaccines
  double _TBVEfficacy;
  //@}
  
  void clearInfection(Infection *iCurrent);
};

#endif
