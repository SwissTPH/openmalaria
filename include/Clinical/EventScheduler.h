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

#ifndef Hmod_ClinicalEventSchduler
#define Hmod_ClinicalEventSchduler

#include "Global.h"
#include "Clinical/ClinicalModel.h"
#include "Clinical/ESCaseManagement.h"
#include "util/AgeGroupInterpolation.h"

#include <boost/unordered_map.hpp>
#include <list>

namespace scnXml {
class HSEventScheduler;
}

namespace OM {
namespace Clinical {

using util::AgeGroupInterpolation;

/** Tracks clinical status (sickness), triggers case management for new events,
 * medicates treatment, determines patient recovery, death and sequelae.
 *
 * Note: there are several variables that only need to be used during a
 * bout. It's possible that memory usage could be reduced by storing them
 * externally in a temporary object during episodes (but unlikely worth doing).
 */
class ClinicalEventScheduler : public ClinicalModel
{
public:
    static void init ();
    static void setParameters (const scnXml::HSEventScheduler& esData);
    static void cleanup ();

    ClinicalEventScheduler (double cF, double tSF);
    ~ClinicalEventScheduler ();
    
    virtual bool notAtRisk();

    virtual void massDrugAdministration(Human& human);

protected:
    virtual void doClinicalUpdate (Human& human, double ageYears);

    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);

private:
    /// Maximum number of timesteps (including first of case) an individual will
    /// remember they are sick before resetting.
    static TimeStep maxUCSeekingMemory;
    /// Length of an uncomplicated case
    static TimeStep uncomplicatedCaseDuration;
    /// Length of a complicated case
    static TimeStep complicatedCaseDuration;
    /// Time-span for which individual is at risk of death in complicated case
    /// minus length of complicated case (must be <= 0)
    static TimeStep extraDaysAtRisk;
    /// First value is probability immediate treatment, second is first +
    /// probability 1-day delay to treatment seeking, etc. Last value must be 1.
    static vector<double> cumDailyPrImmUCTS;

    /// Parameter of S(t) for t > 0
    static double neg_v;
    
    /** Weight model. Currently looks up a weight dependant on age from a table
     * in an entirely deterministic way.
     *
     * @param ageGroupData Age group for weight data
     * @param ageYears Age in years
     * @returns Mass in kg */
    inline double ageToWeight (double ageYears) {
        return weight->eval( ageYears ) * hetWeightMultiplier;
    }
    
    static double hetWeightMultStdDev;
    static double minHetWeightMult;
    static AgeGroupInterpolation* weight;
    
    /** Base log odds of treatment of non-malarial fevers in absense of a
     * malaria diagnostic and irrespective of whether treatment is needed.
     * 
     * In our model, this is logit(P₀), not β₀. */
    static double logOddsAbBase;
    /** Added to log odds treatment when a malaria diagnostic indicates no
     * parasites. Symbol in model: β₁. */
    static double logOddsAbNegTest;
    /** Added to log odds treatment when a malaria diagnostic indicates
     * parasites. Symbol in model: β₂. */
    static double logOddsAbPosTest;
    /** Added to log odds treatment when NMF is categorized as an illness
     * potentially leading to death (Pathogenesis::NEED_ANTIBIOTIC).
     * Symbol in model: β₃. */
    static double logOddsAbNeed;
    /** Added to log odds treatment when given by an informal provider.
     * Symbol in model: β₄. */
    static double logOddsAbInformal;
    /** One minus the efficacy of antibiotic/NMF treatment (i.e. a multiplier
     * for fatality-rate given that the case is treated). */
    static double oneMinusEfficacyAb;
    /** Case fatality rate of non-malaria fevers requiring treatment given that
     * the case is not treated. */
    static AgeGroupInterpolation* severeNmfMortality;
    
    // Note on memory usage: Pathogenesis::State is and enum (an int), so we
    // have a vtable followed by 3 ints, a double and a list. Alignment probably
    // wastes some space.
    /// Current state of sickness
    Pathogenesis::State pgState;

    /** Set to when a bout should start. If TimeStep::simulation equals this, a bout
     * is started (UC & severe behaviour different).
     *
     * Note: medications are not delayed by this. */
    TimeStep caseStartTime;

    /** The individual recovers when TimeStep::simulation >= timeOfRecovery,
     * assuming they didn't die. */
    TimeStep timeOfRecovery;

    /// Time at which last treatment was recieved (for second-case considerations).
    TimeStep timeLastTreatment;

    /// Total parasite density at previous timestep (used during a bout).
    double previousDensity;
    
    /// Multiplies the mean weight for age.
    /// Within PkPd class simply because it's not used elsewhere.
    double hetWeightMultiplier;
    
    /// All pending medications
    list<MedicateData> medicateQueue;
};

}
}
#endif
