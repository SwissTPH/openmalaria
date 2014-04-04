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
      m_inAnyCohort & stream;
      m_subPopExp & stream;
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
  /** Add the human to an intervention component's sub-population for the given
   * duration. A duration of zero implies no effect. */
  void reportDeployment( interventions::ComponentId id, TimeStep duration );
  
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
  
  /** Return true if human is a member of the sub-population.
   * 
   * This is only for use during intervention deployment (see comment on
   * m_subPopExp). */
  inline bool isInSubPop( interventions::ComponentId id )const{
      map<interventions::ComponentId,TimeStep>::const_iterator it = m_subPopExp.find( id );
      if( it == m_subPopExp.end() ) return false;       // no history of membership
      else return it->second > TimeStep::simulation;   // added: has expired?
  }
  /** Return true if human is a member of any cohort. */
  //TODO(monitoring): outputs per cohort, not simply any cohort or everyone
  inline bool isInAnyCohort()const{ return m_inAnyCohort; }
  
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
  
  /** Act on "remove from sub-population on first ..." events.
   *
   * This is only for use during a human update. */
  void removeFirstEvent( interventions::SubPopRemove::RemoveAtCode code );
  
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
  
  ///@brief Cached values used by monitoring
  //@{
  /// Made persistant to save a lookup each timestep (has a significant impact)
  Monitoring::AgeGroup monitoringAgeGroup;
  /// Cache, updated when human is added to or removed from a sub-population
  bool m_inAnyCohort;
  //@}
  
  /// Vaccines
  //TODO: could move TBV code to WHFalciparum, where the efficacy is now used
  interventions::PerHumanVaccine _vaccine;
  
  //!Date of birth, time step since start of warmup
  TimeStep _dateOfBirth;
  
  /// The next continuous distribution in the series
  uint32_t nextCtsDist;
  
  
  /** This lists sub-populations of which the human is a member together with
   * expiry time.
   * 
   * Definition: a human is in a sub-population if that sub-population is
   * listed here, and, at time of intervention deployment the expiry time given
   * here is greater than the current timestep, or during human update the
   * expiry time given here is greater than or equal to the current timestep.
   * NOTE: this discrepancy is because intervention deployment effectively
   * happens at the end of a timestep and we want a duration of 1 timestep to
   * mean 1 intervention deployment (that where the human becomes a member) and
   * 1 human update (the next). */
  map<interventions::ComponentId,TimeStep> m_subPopExp;
};

} }
#endif
