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

#ifndef Hmod_ClinicalModel
#define Hmod_ClinicalModel

#include "Host/Human.h"
#include "Episode.h"
#include <memory>

namespace scnXml{
    class Model;
    class HealthSystem;
}
namespace OM { namespace Clinical {
    using Host::Human;

/** The clinical model models the effects of sickness dependant on malarial
 * parasite densities and administers anti-malaria treatments via the drug
 * model (or in a simpler case, directly clearing infections).
 * 
 * So far, sickness types include uncomplicated and severe malaria cases and
 * non-malaria sickness.
 * 
 * Patient outcomes include full recovery, recovery with sequelae and death.
 * 
 * Reporting includes patient outcome and potentially drug usage and use of
 * RDTs (Rapid Diagnostic Tests) for costing purposes. */\
class ClinicalModel
{
public:
  /// @brief Static functions
  //@{
  /// Initialise whichever model is in use.
  static void init ( const Parameters& parameters, const scnXml::Model& model );
  
    /** Set data for a new health system.
     * 
     */
    static void changeHS (const scnXml::HealthSystem& healthSystem);
  
  /** Return a new ClinicalModel.
   *
   * @param tSF	treatment seeking factor, passed to CaseManagementModel */
  static ClinicalModel* createClinicalModel (double tSF);
  //@}
  
  /// Destructor
  virtual ~ClinicalModel ();
  
  /** Kills the human if ageTimeSteps reaches the simulation age limit.
   *
   * @returns True if the human has been killed by some means. The clinical
   *	model now tracks this status. */
  bool isDead (TimeStep ageTimeSteps);
  
  /** Run main part of the model: determine the sickness status and any
   * treatment for the human.
   * 
   * @param ageYears Age of human.
   * @param ageTimeSteps Age of human (used to test if 1 timestep old) */
  void update (Human& human, double ageYears, TimeStep ageTimeSteps);
  
  /** For infants, updates the infantIntervalsAtRisk and potentially
   * infantDeaths arrays. */
  void updateInfantDeaths (TimeStep ageTimeSteps);
  
  /** Used with IPT within host model to potentially avoid further reports:
   * The four timesteps after a bout are not at risk of a further bout since
   * if one occured it would be considered the same bout.
   *
   * Only supported by immediate outcomes model. */
  virtual bool notAtRisk() =0;
  
    virtual void massDrugAdministration( Human& human,
        Monitoring::ReportMeasureI screeningReport,
        Monitoring::ReportMeasureI drugReport ) =0;
  
  /// Force all pending summaries to be reported. Should only be called when
  /// class is about to be destroyed anyway to avoid affecting output.
  inline void flushReports (){
      latestReport.flush();
  }
  
  /// Checkpointing
  template<class S>
  void operator& (S& stream) {
      checkpoint (stream);
  }
  
protected:
  /// Constructor.
  ClinicalModel ();
  
  /** Update for clinical model - new pathogenesis status, treatment, etc.
   *
   * @param withinHostModel WithinHostModel of human.
   * @param hostTransmission per-host transmission data of human.
   * @param ageYears Age of human.
   * @param ageGroup Survey age group of human. */
  virtual void doClinicalUpdate (Human& human, double ageYears) =0;
  
  virtual void checkpoint (istream& stream);
  virtual void checkpoint (ostream& stream);
  
  /** Last episode; report to survey pending a new episode or human's death. */
  Episode latestReport;
  
  /** @brief Positive values of _doomed variable (codes). */
  enum {
    DOOMED_TOO_OLD = 1,		///< died because reached age limit
    DOOMED_COMPLICATED = 4,	///< died from severe malaria or malaria with a coinfection
    DOOMED_NEONATAL = 6,	///< died due to mother's malaria infection
    DOOMED_INDIRECT = 7,	///< died indirectly from malaria (after a delay)
  };
  
  /** Can indicate that the individual is dead or about to die.
   *
   * If _doomed < 0, the individual is doomed to die.
   * 
   * If _doomed > 0, the individual is dead, and will be removed from the
   * population at the beginning of the next timestep. NOTE: why not
   * immediately? See above enum for positive values used. */
  int _doomed;
};

} }
#endif
