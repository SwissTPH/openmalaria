/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2014 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2014 Liverpool School Of Tropical Medicine
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
#include "Transmission/PerHost.h"
#include "InfectionIncidenceModel.h"
#include "Monitoring/Survey.h"
#include "Monitoring/AgeGroup.h"
#include "interventions/HumanComponents.h"
#include "util/checkpoint_containers.h"
#include <set>

namespace scnXml {
    class Scenario;
}
namespace OM {
namespace Transmission {
    class TransmissionModel;
}
namespace Clinical {
    class ClinicalModel;
}
namespace WithinHost {
    class WHInterface;
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
      _dateOfBirth & stream;
      _vaccine & stream;
      nextCtsDist & stream;
      cohorts & stream;
      lastDeployments & stream;
  }
  //@}
  
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
  /** Mark a certain intervention component as being deployed now. */
  inline void updateLastDeployed( interventions::ComponentId id ){
      lastDeployments[id] = TimeStep::simulation;
  }
  
  /** Determines for the purposes of cumulative deployment whether a component is
   * still current.
   * 
   * @returns true if the intervention component should be re-deployed (too old)
   */
  bool needsRedeployment( interventions::ComponentId cumCuvId, TimeStep maxAge );
  
  /// Resets immunity
  void clearImmunity();
  
  /// Infect the human (with an imported infection).
  void addInfection();
  //@}
  
  /// @brief Small functions
  //@{
    /** Get the age in time steps, based on current TimeStep::simulation. */
    inline TimeStep getAgeInTimeSteps() const{
        return TimeStep::simulation - _dateOfBirth;
    }
    /** Get the age in years, based on current TimeStep::simulation. */
    inline double getAgeInYears() const{
        return (TimeStep::simulation - _dateOfBirth).inYears();
    }
  
  //! Returns the date of birth
  inline TimeStep getDateOfBirth() {return _dateOfBirth;}
  
  /** Return true if human is a member of the cohort.
   * 
   * Note: the maximum value of size_t is considered to be the population, of
   * which every human is a member. */
  inline bool isInCohort( interventions::ComponentId id )const{
      return id == interventions::ComponentId_pop || cohorts.count( id ) > 0;
  }
  /** Return true if human is a member of any cohort. */
  //TODO(monitoring): outputs per cohort, not simply any cohort or everyone
  inline bool isInAnyCohort()const{ return cohorts.size() > 0; }
  
  /// Return the index of next continuous intervention to be deployed
  inline uint32_t getNextCtsDist()const{ return nextCtsDist; }
  /// Increment then return index of next continuous intervention to deploy
  inline uint32_t incrNextCtsDist() {
      nextCtsDist += 1;
      return nextCtsDist;
  }
  //@}
  
  //! Summarize the state of a human individual.
  void summarize();
  
  /** Add human to a cohort (assumes cohort mode is active).
   *
   * Also makes sure inter-survey stats will only be
   * summed from this point onwards (i.e. removes data accumulated between
   * last time human was reported or birth and now). */
  void addToCohort ( interventions::ComponentId, Monitoring::ReportMeasureI reportMeasure );
  
  /** Remove from a cohort. As with addToCohort, deals with reporting.
   *
   * Can be safely called when human is not in cohort. */
  void removeFromCohort( interventions::ComponentId );
  
  /** Act on remove-from-cohort-on-first-xyz events. */
  void removeFromCohorts( interventions::Cohort::RemoveAtCode code );
  
  /// Flush any information pending reporting. Should only be called at destruction.
  void flushReports ();
  
  ///@brief Access to sub-models
  //@{
  /// The WithinHostModel models parasite density and immunity
  inline const WithinHost::WHInterface& getWithinHostModel () const{
      return *withinHostModel;
  }
  
  /// Get monitoring age group
  inline const Monitoring::AgeGroup& getMonitoringAgeGroup() const{
      return monitoringAgeGroup;
  }
  
  inline interventions::PerHumanVaccine& getVaccine(){ return _vaccine; }
  inline const interventions::PerHumanVaccine& getVaccine() const{ return _vaccine; }
  
  inline Clinical::ClinicalModel& getClinicalModel() {
      return *clinicalModel;
  }
  //@}
  
  
  ///@name static public
  //@{
  static void initHumanParameters (const OM::Parameters& parameters, const scnXml::Scenario& scenario);
  
  static void clear();
  //@}
    
public:
  /** @brief Models
   *
   * These contain various sub-models used by Humans. */
  //@{
  /// Contains per-species vector data (VectorModel only).
  Transmission::PerHost perHostTransmission;
  
  /// The WithinHostModel models parasite density and immunity
  WithinHost::WHInterface *withinHostModel;
  
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
  //TODO: could move TBV code to WHFalciparum, where the efficacy is now used
  interventions::PerHumanVaccine _vaccine;
  
  //!Date of birth, time step since start of warmup
  TimeStep _dateOfBirth;
  
  /// The next continuous distribution in the series
  uint32_t nextCtsDist;
  
  
  //TODO(performance): are dynamic maps/sets the best approach or is using vector or boost::dynamic_bitset better?
  /// The set of cohorts (intervention indexes) to which this human is a member
  set<interventions::ComponentId> cohorts;
  /// Last deployment times of intervention components by component id
  map<interventions::ComponentId,TimeStep> lastDeployments;
};

} }
#endif
