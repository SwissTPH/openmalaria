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
#include "Global.h"
#include "Host/Vaccine.h"
#include "Host/ContinuousIntervention.h"
#include "Transmission/PerHostTransmission.h"
#include "InfectionIncidenceModel.h"
#include "WithinHost/WithinHostModel.h"
#include "Monitoring/Survey.h"

namespace OM {
    namespace Transmission {
	class TransmissionModel;
    }
    namespace Clinical {
	class ClinicalModel;
    }
namespace Host {

/** Interface to all sub-models storing data per-human individual.
 *
 * Still contains some data, but most is now contained in sub-models. */
class Human {
public:
  
  /// @brief Construction and destruction, checkpointing
  //@{
  /** Initialise all variables of a human datatype.
   * 
   * \param tm Transmission model reference (to initialize TM code)
   * \param dateOfBirth date of birth in time steps (equal to simulationTime,
   *	except for initial population set up)
   * \param simulationTime Simulation timestep */
  Human(Transmission::TransmissionModel& tm, int dateOfBirth, int simulationTime);

  /** Destructor
   * 
   * Note: this destructor does nothing in order to allow shallow copying of a
   * Human into the population list. Human::destroy() does the real freeing and
   * must be called explicitly. */
  ~Human() {}
  
  /// The real destructor
  void destroy();
  
  /// Checkpointing
  template<class S>
  void operator& (S& stream) {
      perHostTransmission & stream;
      // In this case these pointers each refer to one element not stored/pointed
      // from elsewhere, so this checkpointing technique works.
      (*infIncidence) & stream;
      (*withinHostModel) & stream;
      (*clinicalModel) & stream;
      monitoringAgeGroup & stream;
      ageGroupData & stream;
      ctsIntervention & stream;
      _ylag & stream;
      _dateOfBirth & stream;
      _vaccine & stream;
      _probTransmissionToMosquito & stream;
      _inCohort & stream;
  }
  //@}
  
  /** @brief Per timestep update functions
   * Functions called by Population, Simulation or Summary */
  //@{
  /// Update ageGroupData. Should happen before anything using this on  the timestep.
  inline void updateAgeGroupData() {
      ageGroupData.update( getAgeInYears() );
  }
      
  /** Main human update.
   *
   * @param simulationTime Time-step
   * @param transmissionModel Pointer to transmission data.
   * @param doUpdate If false, returns immediately after is-dead check.
   * @returns True if the individual is dead (too old or otherwise killed).
   */
  bool update(int simulationTime, Transmission::TransmissionModel* transmissionModel, bool doUpdate);
  //@}
  
  ///@brief Deploy "intervention" functions
  //@{
  /// Asks the clinical model to deal with this
  void massDrugAdministration ();
  
  /// Vaccinate & report mass vaccination
  void massVaccinate ();
  /// If individual hasn't dropped out, vaccinate & report EPI
  void ctsVaccinate ();
  
  void IPTiTreatment ();
  void deployIptDose ();
  
  /// Give human a new ITN via mass deployment
  void massITN ();
  /// Give a human a new ITN through EPI
  void ctsITN ();
  /// Give human a new IRS through mass deployment
  void massIRS ();
  /// Give human a new VA intervention through mass deployment
  void massVA ();
  
  /// Resets immunity
  inline void immuneSuppression() {
      withinHostModel->immuneSuppression();
  }
  
  /// Infect the human (with an imported infection).
  void addInfection();
  
  /// Add PEV and remove TBV (vaccines) from human
  inline void R_0Vaccines() {
      _vaccine.setPEV( 1.0 );
      _vaccine.setTBV( 0.0 );
  }
  //@}
  
  /// @brief Small functions
  //@{
  /// Return ageGroupData
  inline const AgeGroupData getAgeGroupData() const{
      return ageGroupData;
  }
      
  /// Get the survey age-group. Constant-time; returns result of last update.
  inline const Monitoring::AgeGroup ageGroup() const {
      return monitoringAgeGroup;
  }
  
  //! Get the age in years, based on current simulationTime.
  double getAgeInYears() const;
  
  //! Returns the date of birth
  inline int getDateOfBirth() {return _dateOfBirth;}
  
  /** Does the Human have a detectible infection? */
  inline bool detectibleInfection () const {
    return withinHostModel->parasiteDensityDetectible();
  }
  
  inline bool getInCohort(){ return _inCohort; }
  //@}
  
  //! Summarize the state of a human individual.
  void summarize();
  
  /** Add human to a cohort (assumes cohort mode is active).
   *
   * Also makes sure inter-survey stats will only be
   * summed from this point onwards (i.e. removes data accumulated between
   * last time human was reported or birth and now). */
  void addToCohort ();
  
  /// Flush any information pending reporting. Should only be called at destruction.
  void flushReports ();
  
  /// Calculate chance of a biting mosquito becoming infected
  //TODO: per genotype? (for Tiago's resistance modelling)
  inline double probTransmissionToMosquito() const {
    return _probTransmissionToMosquito;
  }
  
  ///@brief Access to sub-models
  //@{
  /// The WithinHostModel models parasite density and immunity
  inline const WithinHost::WithinHostModel& getWithinHostModel (){
      return *withinHostModel;
  }
  //@}
  
  
  ///@name static public
  //@{
  static void initHumanParameters ();
  
  static void clear();
  //@}
  
private:
    void updateInfection(Transmission::TransmissionModel*, double ageYears);
    
    /*! Apply interventions to this human if eligible. Calculate the remaining
    efficacy of the latest vaccination if vaccinated before */
    void updateInterventionStatus();
    
    double calcProbTransmissionToMosquito() const;
    
    void clearInfection(WithinHost::Infection *iCurrent);
    
    
public:
  /** @brief Models
   *
   * These contain various sub-models used by Humans. */
  //@{
  /// Contains per-species vector data (VectorTransmission only).
  Transmission::PerHostTransmission perHostTransmission;
  
  /// The InfectionIncidenceModel translates per-host EIR into new infections
  InfectionIncidenceModel *infIncidence;
  
private:
  /// The WithinHostModel models parasite density and immunity
  WithinHost::WithinHostModel *withinHostModel;
  
  /** The ClinicalModel encapsulates pathogenesis (sickness status),
   * case management (medicating drugs)
   * and clinical outcomes (morbidity, reporting). */
  Clinical::ClinicalModel *clinicalModel;
  //@}
  
  /// Made persistant to save a lookup each timestep (has a significant impact)
  Monitoring::AgeGroup monitoringAgeGroup;
  
  /// Age group for compiled-in data (made persistant as for monitoringAgeGroup).
  AgeGroupData ageGroupData;
  
  /// Continuous intervention deployment
  ContinuousIntervention ctsIntervention;
  
  /// Vaccines
  PerHumanVaccine _vaccine;
  
  /** Total asexual blood stage density over last 20 days (need samples from 10, 15 and 20 days ago)
   *
   * _ylag[simulationTime % _ylagLen] corresponds to current timestep. */
  vector<double> _ylag;
  /// Length of _ylag array. Wouldn't have to be dynamic if Global::interval was known at compile-time.
  /// set by initHumanParameters
  static int _ylagLen;
  
  //!Date of birth, time step since start of warmup
  int _dateOfBirth;
  
  /// True if human is included in a cohort.
  bool _inCohort;
  
  /// Cached value of calcProbTransmissionToMosquito; checkpointed
  double _probTransmissionToMosquito;
};

} }
#endif
