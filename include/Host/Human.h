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
#include "Transmission/PerHostTransmission.h"
#include "InfectionIncidenceModel.h"
#include "WithinHost/WithinHostModel.h"
#include "Monitoring/Surveys.h"

namespace OM {
    namespace Transmission {
	class TransmissionModel;
    }
    namespace Clinical {
	class ClinicalModel;
    }
    class Population;
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
   * \param dateOfBirth date of birth in time steps (equal to TimeStep::simulation,
   *	except for initial population set up)
   */
  Human(Transmission::TransmissionModel& tm, TimeStep dateOfBirth);

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
      _ylag & stream;
      _dateOfBirth & stream;
      _vaccine & stream;
      _probTransmissionToMosquito & stream;
      _inCohort & stream;
      nextCtsDist & stream;
  }
  //@}
      
  /** Main human update.
   *
   * @param transmissionModel Pointer to transmission data.
   * @param doUpdate If false, returns immediately after is-dead check.
   * @returns True if the individual is dead (too old or otherwise killed).
   */
  bool update(const OM::Population&, Transmission::TransmissionModel* transmissionModel, bool doUpdate);
  //@}
  
  ///@brief Deploy "intervention" functions
  //@{
  /// Asks the clinical model to deal with this
  void massDrugAdministration (const OM::Population&);
  
  /// Vaccinate & report mass vaccination
  void massVaccinate (const OM::Population&);
  /// If individual hasn't dropped out, vaccinate & report EPI
  void ctsVaccinate (const OM::Population&);
  
  void IPTiTreatment (const OM::Population&);
  void deployIptDose (const OM::Population&);
  
  /// Give human a new ITN via mass deployment
  void massITN (const OM::Population&);
  /// Give a human a new ITN through EPI
  void ctsITN (const OM::Population&);
  /// Give human a new IRS through mass deployment
  void massIRS (const OM::Population&);
  /// Give human a new VA intervention through mass deployment
  void massVA (const OM::Population&);
  
  /// Resets immunity
  inline void immuneSuppression(const OM::Population&) {
      withinHostModel->immuneSuppression();
  }
  
  /// Infect the human (with an imported infection).
  void addInfection();
  
  /// Add PEV and remove TBV (vaccines) from human
  inline void R_0Vaccines() {
      _vaccine.setInitialPEV( 1.0 );
      _vaccine.setInitialTBV( 0.0 );
  }
  //@}
  
  ///@brief Functions to check coverage by interventions
  //@{
    bool hasVaccineProtection(TimeStep maxInterventionAge) const;
    bool hasIPTiProtection(TimeStep maxInterventionAge) const;
    bool hasITNProtection(TimeStep maxInterventionAge) const;
    bool hasIRSProtection(TimeStep maxInterventionAge) const;
    bool hasVAProtection(TimeStep maxInterventionAge) const;
  //@}
  
  /// @brief Small functions
  //@{
  //! Get the age in years, based on current TimeStep::simulation.
  double getAgeInYears() const;
  
  //! Returns the date of birth
  inline TimeStep getDateOfBirth() {return _dateOfBirth;}
  
  /** Does the Human have a detectible infection? */
  inline bool detectibleInfection () const {
    return withinHostModel->parasiteDensityDetectible();
  }
  
  // crux for timed deployment as intervention up to some limit:
  inline bool getInCohort(TimeStep)const{ return _inCohort; }
  inline bool getInCohort()const{ return _inCohort; }
  //@}
  
  /// Return the current survey to use (depends on survey time and whether or
  /// not individual is in the cohort).
  Monitoring::Survey& getSurvey() const{
      return Monitoring::Surveys.getSurvey( _inCohort );
  }
  
  //! Summarize the state of a human individual.
  void summarize();
  
  /** Add human to a cohort (assumes cohort mode is active).
   *
   * Also makes sure inter-survey stats will only be
   * summed from this point onwards (i.e. removes data accumulated between
   * last time human was reported or birth and now). */
  void addToCohort (const OM::Population&);
  
  /** Remove from cohort. As with addToCohort, deals with reporting.
   *
   * Can be safely called when human is not in cohort. */
  void removeFromCohort();
  
  /// Flush any information pending reporting. Should only be called at destruction.
  void flushReports ();
  
  /** Calculate chance of a biting mosquito becoming infected
   * 
   * Before update() is called, this refers to the value for this time-step.
   * Afterwards, it refers to the value for the *next* time-step. */
  //TODO: per genotype? (for LSTM's spread of resistance modelling)
  inline double probTransmissionToMosquito() const {
    return _probTransmissionToMosquito;
  }
  
  ///@brief Access to sub-models
  //@{
  /// The WithinHostModel models parasite density and immunity
  inline const WithinHost::WithinHostModel& getWithinHostModel () const{
      return *withinHostModel;
  }
  
  /// Get monitoring age group
  inline const Monitoring::AgeGroup& getMonitoringAgeGroup() const{
      return monitoringAgeGroup;
  }
  
  inline const PerHumanVaccine& getVaccine() const{
      return _vaccine;
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
    void updateInterventionStatus(const OM::Population&);
    
    double calcProbTransmissionToMosquito() const;
    
    void clearInfection(WithinHost::Infection *iCurrent);
    
    
public:
  /** @brief Models
   *
   * These contain various sub-models used by Humans. */
  //@{
  /// Contains per-species vector data (VectorTransmission only).
  Transmission::PerHostTransmission perHostTransmission;
  
  /// The WithinHostModel models parasite density and immunity
  WithinHost::WithinHostModel *withinHostModel;
  
private:
  /// The InfectionIncidenceModel translates per-host EIR into new infections
  InfectionIncidenceModel *infIncidence;
  
  /** The ClinicalModel encapsulates pathogenesis (sickness status),
   * case management (medicating drugs)
   * and clinical outcomes (morbidity, reporting). */
  Clinical::ClinicalModel *clinicalModel;
  //@}
  
  /// Made persistant to save a lookup each timestep (has a significant impact)
  Monitoring::AgeGroup monitoringAgeGroup;
  
  /// Vaccines
  PerHumanVaccine _vaccine;
  
  /** Total asexual blood stage density over last 20 days (need samples from 10, 15 and 20 days ago)
   *
   * _ylag[TimeStep::simulation % _ylagLen] corresponds to current time step. */
  vector<double> _ylag;
  /// Length of _ylag array. Wouldn't have to be dynamic if Global::interval was known at compile-time.
  /// set by initHumanParameters
  static int _ylagLen;
  
  //!Date of birth, time step since start of warmup
  TimeStep _dateOfBirth;
  
  /// The next continuous distribution in the series
  uint32_t nextCtsDist;
  
  /// True if human is included in a cohort.
  bool _inCohort;
  
  /// Cached value of calcProbTransmissionToMosquito; checkpointed.
  ///
  /// After human update, contains value corresponding to next time-step.
  double _probTransmissionToMosquito;
  
public: //lazy: give read access to these
  /// Remove from cohort as soon as individual has patent parasites?
  static bool cohortFirstInfectionOnly;
  /// Remove from cohort as soon as individual receives treatment?
  static bool cohortFirstTreatmentOnly;
  /// Remove from cohort as soon as individual gets sick (any sickness)?
  static bool cohortFirstBoutOnly;
};

} }
#endif
