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

#ifndef Hmod_ClinicalImmediateOutcomes
#define Hmod_ClinicalImmediateOutcomes

#include "WithinHost/Pathogenesis/State.h"
#include "Clinical/ClinicalModel.h"
#include "WithinHost/WHInterface.h"

namespace OM {
namespace Clinical {

namespace Regimen {
/** Regimen: UC / UC2 / SEVERE.
 *
 * Note: values used in array lookups, so are important. */
enum Type {
    UC = 0,         // first line
    UC2 = 1,                // second line
    SEVERE = 2,     // third line
    NUM = 3,
};
}

/** This implementation of the model is intended to use the old case-management
 * model with immediate outcomes of clinical events (immediate recovery with
 * total parasite clearance or immediate death). */
class ClinicalImmediateOutcomes : public ClinicalModel
{
public:
    /** Initialises parameters, loading from XML data. */
    static void initParameters ();

    /** Load health system data from initial data or an intervention's data (both from XML).
     * (Re)loads all data affected by this healthSystem element. */
    static void setHealthSystem (const scnXml::HealthSystem& healthSystem);

    ClinicalImmediateOutcomes (double tSF);
    ~ClinicalImmediateOutcomes ();

    virtual bool notAtRisk() {
        return (TimeStep::simulation-_tLastTreatment >= TimeStep(1) &&
                TimeStep::simulation-_tLastTreatment <= TimeStep(4));
    }

    virtual void massDrugAdministration( Human& human,
        Monitoring::ReportMeasureI screeningReport,
        Monitoring::ReportMeasureI drugReport );

protected:
    virtual void doClinicalUpdate (Human& human, double ageYears);

    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);

private:
    /** Called when a non-severe/complicated malaria sickness occurs. */
    void uncomplicatedEvent(Human& human, Episode::State pgState);

    /** Called when a severe/complicated (with co-infection) malaria sickness occurs.
     *
     * Note: sets doomed = 4 if patient dies. */
    void severeMalaria(Human& human, Episode::State pgState, double ageYears, int& doomed);

    /** Timestep of the last treatment (TIMESTEP_NEVER if never treated). */
    TimeStep _tLastTreatment;

    //! treatment seeking for heterogeneity
    double _treatmentSeekingFactor;
    
    /// Calculate _probGetsTreatment, _probParasitesCleared and _cureRate.
    static void setParasiteCaseParameters (const scnXml::HSImmediateOutcomes& healthSystem);

    //BEGIN Static parameters, set by setHealthSystem()
    // These parameters are reset via a setHealthSystem call on checkpoint
    // load rather than checkpointed.

    static double probGetsTreatment[Regimen::NUM];
    static double probParasitesCleared[Regimen::NUM];
    static double cureRate[Regimen::NUM];
    static WithinHost::TreatmentId treatments[Regimen::NUM];
    //END
};

}
}
#endif
