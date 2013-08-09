/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2013 Swiss Tropical and Public Health Institute 
 * Copyright (C) 2005-2013 Liverpool School Of Tropical Medicine
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
#include "Transmission/PerHost.h"
#include "InfectionIncidenceModel.h"
#include "WithinHost/WithinHostModel.h"
#include "Monitoring/Surveys.h"
#include "interventions/Interventions.h"
#include <set>

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
      nextCtsDist & stream;
      cohorts & stream;
      lastDeployments & stream;
  }
  //@}
  
  /** Update infectiousness to mosquitoes. */
  void updateInfectiousness();
  
  /** Main human update.
   *
   * @param transmissionModel Pointer to transmission data.
   * @param doUpdate If false, returns immediately after is-dead check.
   * @returns True if the individual is dead (too old or otherwise killed).
   */
  bool update(OM::Transmission::TransmissionModel* transmissionModel, bool doUpdate);
  //@}
  
  ///@brief Deploy "intervention" functions
  //@{
  /** Deploys an intervention.
   * 
   * The purpose of this function is principally to allow recording of when
   * each intervention effect is deployed (for cumulative deployment).
   * 
   * The number of wrapper functions in intervention deployment could probably
   * be reduced somehow (TODO). */
  void deploy( const interventions::HumanInterventionEffect& effect, interventions::Deployment::Method method );
  
  /** Determines for the purposes of cumulative deployment whether an effect is
   * still current.
   * 
   * @returns true if the intervention effect should be re-deployed (too old)
   */
  bool needsRedeployment( size_t effect_index, TimeStep maxAge );
  
  /// Mass/EPI vaccination (note: extra checks may still prevent EPI vaccination)
  void deployVaccine( interventions::Deployment::Method method, Vaccine::Types type );
  /// Mass/continuous IPT deployment (note: these don't have exactly the same effect)
  void deployIPT( interventions::Deployment::Method method );
  /// Mass/continuous ITN deployment
  void deployITN( interventions::Deployment::Method method, Transmission::TransmissionModel& transmissionModel );
  /// Mass/continuous IRS deployment
  void deployIRS( interventions::Deployment::Method method, Transmission::TransmissionModel& transmissionModel );
  
  /// Resets immunity
  inline void clearImmunity() {
      withinHostModel->clearImmunity();
  }
  
  /// Infect the human (with an imported infection).
  void addInfection();
  
  /// Add PEV and remove TBV (vaccines) from human
  inline void R_0Vaccines() { _vaccine.specialR_0(); }
  
  /** Report deployment of an intervention to this human. */
  void reportDeployment( interventions::Effect::Type type, interventions::Deployment::Method method ) const;
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
  
  /** Return true if human is a member of the cohort.
   * 
   * Note: the maximum value of size_t is considered to be the population, of
   * which every human is a member. */
  inline bool isInCohort( size_t index )const{
      return index == numeric_limits<size_t>::max() || cohorts.count( index ) > 0;
  }
  /** Return true if human is a member of any cohort.
   *
   * (TODO: code should be updated not to use this) */
  inline bool isInAnyCohort()const{ return cohorts.size() > 0; }
  
  /// Return the index of next continuous intervention to be deployed
  inline uint32_t getNextCtsDist()const{ return nextCtsDist; }
  /// Increment then return index of next continuous intervention to deploy
  inline uint32_t incrNextCtsDist() {
      nextCtsDist += 1;
      return nextCtsDist;
  }
  //@}
  
  /// Return the current survey to use (depends on survey time and whether or
  /// not individual is in the cohort).
  Monitoring::Survey& getSurvey() const{
      //TODO: inappropriate use of cohorts
      return Monitoring::Surveys.getSurvey( cohorts.size() > 0 );
  }
  
  //! Summarize the state of a human individual.
  void summarize();
  
  /** Add human to a cohort (assumes cohort mode is active).
   *
   * Also makes sure inter-survey stats will only be
   * summed from this point onwards (i.e. removes data accumulated between
   * last time human was reported or birth and now). */
  void addToCohort ( size_t index );
  
  /** Remove from a cohort. As with addToCohort, deals with reporting.
   *
   * Can be safely called when human is not in cohort. */
  void removeFromCohort( size_t index );
  /** As removeFromCohort, but for all cohorts. */ //TODO: check usages of this
  void removeFromAllCohorts();
  
  /// Flush any information pending reporting. Should only be called at destruction.
  void flushReports ();
  
  /** Return the infectiousness of this human to biting mosquitoes.
   * 
   * Returns the value for the last time-step on which updateInfectiousness
   * has been called. */
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
  
  inline const Clinical::ClinicalModel& getClinicalModel() const{
      return *clinicalModel;
  }
  //@}
  
  
  ///@name static public
  //@{
  static void initHumanParameters ();
  
  static void clear();
  //@}
  
private:
    void updateInfection(Transmission::TransmissionModel*, double ageYears);
    
    void clearInfection(WithinHost::Infection *iCurrent);
    
    
public:
  /** @brief Models
   *
   * These contain various sub-models used by Humans. */
  //@{
  /// Contains per-species vector data (VectorModel only).
  Transmission::PerHost perHostTransmission;
  
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
  
  /** Total asexual blood stage density over last 20 days (uses samples from
   * 10, 15 and 20 days ago).
   *
   * _ylag[mod(TimeStep::simulation, _ylagLen)] corresponds to the density from the
   * previous time step (once updateInfection has been called). */
  vector<double> _ylag;
  /// Length of _ylag array. Wouldn't have to be dynamic if Global::interval was known at compile-time.
  /// set by initHumanParameters
  static int _ylagLen;
  
  //!Date of birth, time step since start of warmup
  TimeStep _dateOfBirth;
  
  /// The next continuous distribution in the series
  uint32_t nextCtsDist;
  
  /// Cached value of calcProbTransmissionToMosquito; checkpointed
  double _probTransmissionToMosquito;
  
  //TODO: are dynamic maps/sets the best approach or is using vector or boost::dynamic_bitset better?
  /// The set of cohorts (intervention indexes) to which this human is a member
  set<size_t> cohorts;
  /// Last deployment times of intervention effects by effect index
  map<size_t,TimeStep> lastDeployments;
  
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
