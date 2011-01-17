/*

  This file is part of OpenMalaria.
 
  Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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

#ifndef Hmod_LSTMDrugType
#define Hmod_LSTMDrugType

#include <string>
#include <deque>
#include <map>
#include <vector>
#include <set>
#include <cassert>
#include "Global.h"

using namespace std;

namespace scnXml{
    class DrugDescription;
    class Drug;
    class Allele;
}
namespace OM { namespace PkPd {

class LSTMDrugType;

/** Per drug, per allele parameters and functions to calculate drug factors and
 * concentrations. */
class LSTMDrugAllele {
    struct Cache {
        Cache( double c, double d, double r );
        
        // hash, used as a way of organizing in a map
        size_t hash;
        // inputs:
        double C0, duration, rate;
        // cached outputs:
        double C1, drugFactor;
        
        // check inputs are equal (used to assert relation):
        bool operator== (const Cache& rhs) const{
            return C0 == rhs.C0 && duration == rhs.duration && rate == rhs.rate;
        }
    };
    struct CacheCompare {
        bool operator()(const Cache& x,const Cache& y) const {
            return x.hash>y.hash;
        }
    };
    typedef std::set<Cache,CacheCompare> CachedIV;
    mutable CachedIV cachedIV;
    
    /// Slope of the dose response curve (no unit)
    double slope;
    /// Maximal drug killing rate per day / (elimination_rate_constant * slope) (no unit)
    double power;
    /// Concentration with 50% of the maximal parasite killing to-the-power-of slope ((mg/l)^slope)
    double IC50_pow_slope;
    /// Maximal drug killing rate per day
    double max_killing_rate;    //TODO: keep this and power?
    
public:
    LSTMDrugAllele( const scnXml::Allele& allele, double elimination_rate_constant );
    
    /** Calculate a survival factor induced by a drug already in the blood.
     * It is expected that no drug doses are taken over the period for which
     * this function calculates a drug factor.
     * 
     * @param drug Reference to per-drug data
     * @param C0 Concentration of drug in blood at start of period. Will be
     *  updated to correct concentration at end of period.
     * @param duration Length of IV in days.
     */
    double calcFactor( const LSTMDrugType& drug, double& C0, double duration ) const;
    
    /** Calculate a survival factor over the course of an intravenous transfusion.
     * No other drug administration should happen during this time span.
     *
     * @param drug Reference to per-drug data
     * @param C0 Concentration of drug in blood at start of IV. Will be
     *  updated to correct concentration at end of IV.
     * @param duration Length of IV in days.
     * @param rate Rate of drug administration (mg/kg/day)
     */
    double calcFactorIV( const LSTMDrugType& drug, double& C0, double duration, double rate ) const;
};
    
    
/** Information about each (type of) drug (rather than each use of a drug).
 *
 * Static data contains a list of all available drug types.
 * 
 * No DrugType data is checkpointed, because it is loaded by init() from XML
 * data. (Although if it cannot be reproduced by reloading it should be
 * checkpointed.) */
class LSTMDrugType {
public:
    ///@brief Static functions
    //@{
    /** Initialise the drug model. Called at start of simulation. */
    static void init (const scnXml::DrugDescription& data);
    /// Remove set-up drugs. (Must be called before init can be re-called.)
    static void cleanup ();
    
    //! Adds a new drug type to the list
    static void addDrug(const LSTMDrugType* drug);
    
    /** Find a DrugType by its abbreviation, and create a new Drug from that.
     *
     * Throws if the drug isn't found, so you can rely on it returning a valid
     * drug if it returns (doesn't throw). */
    static const LSTMDrugType& getDrug(string abbreviation);
    
    /// Return a new proteome ID
    static uint32_t new_proteome_ID ();
    //@}
    
    
    ///@brief Non static (per instance) functions
    //@{
    /** Create a new DrugType.
     *
     * @param drugData Scenario data for this drug (PK params, PD params per allele)
     * @param bit_start Next bit of infection's proteome_id available (see allele_rshift).
     */
    LSTMDrugType (const scnXml::Drug& drugData, uint32_t& bit_start);
    ~LSTMDrugType ();
    
    inline const string& getAbbreviation() const{
        return abbreviation;
    }
    inline double getVolumeOfDistribution() const{
        return vol_dist;
    }
    inline double getNegligibleConcentration() const{
        return negligible_concentration;
    }
    
    /** Return reference to correct drug-allele data. */
    const LSTMDrugAllele& getAllele( uint32_t proteome_ID ) const;
    
    /** Decay concentration C0 over time duration (days) assuming no
     * administration during this time. */
    void updateConcentration( double& C0, double duration ) const;
    /** Update concentration C0 over time duration (days) assuming an
     * intravenous infusion at _rate rate (mg/kg/day) and no
     * administration during this time. */
    void updateConcentrationIV( double& C0, double duration, double rate ) const;
    //@}
  
private:
    // non-copyable (due to allocation of members of drugAllele)
    LSTMDrugType( const LSTMDrugType& );
    
    // The list of available drugs. Not checkpointed; should be set up by init().
    typedef map<const string,const LSTMDrugType*> Available;
    static Available available;
    
    //! The drug abbreviated name, used for registry lookups.
    string abbreviation;
    
    /** Allele information is stored as a uint32_t in infection. Denote this p_id,
     * then we use ((p_id >> allele_rshift) & allele_mask) as an index in
     * PD_params for the allele.
     * 
     * This does restrict the number of alleles, for all drugs, that can be
     * represented, so might have to be changed or extended. */
    uint32_t allele_rshift, allele_mask;
    
    /*PD parameters required - varies with infection genotype*/
    vector<LSTMDrugAllele*> drugAllele;
    
    /*PK parameters required - varies with humans age and severity of disease*/
    /** Concentration, below which drug is deemed not to have an effect and is
     * removed for performance reasons. (mg/l) */
    double negligible_concentration;
    /** Terminal elimination rate constant (negated). Found using
     * ln(2)/half_life. (1 / days) */
    double neg_elimination_rate_constant;
    /// Volume of distribution (l/kg)
    double vol_dist;
    
    /* Resistance data */
    /** Cumulative initial frequencies of each allele. Length and indicies
     * correspond to drugAllele vector.
     * 
     * Independant of frequencies of alleles at other loci (for other drugs). */
    vector<double> cumInitialFreq;
    
    // Allow LSTMDrug to access private members
    friend class LSTMDrugAllele;
    friend inline double drugEffect (const LSTMDrugType& drugType, double& concentration, double duration, double weight_kg, double dose_mg);
};

} }
#endif