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

#ifndef Hmod_LSTMTreatments
#define Hmod_LSTMTreatments

#include <string>
#include <deque>
#include <vector>
#include <set>
#include <cassert>
#include <memory>
#include <map>
#include <boost/unordered_set.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include "Global.h"

using namespace std;

namespace scnXml{
    class Treatments;
}
namespace OM { namespace PkPd {
    
/** Information about treatment courses: schedules and dosage tables.
 * 
 * A schedule is a list of times to administrer a drug, along with base drug
 */
class LSTMDrugType {
public:
    ///@brief Static functions
    //@{
    /** Initialise the drug model. Called at start of simulation. */
    static void init (const scnXml::Pharmacology& data);
    /// Remove set-up drugs. (Must be called before init can be re-called.)
    static void cleanup ();
    
    /** Adds a new drug type to the list. This function becomes responsible
     * for deleting the object on exit. */
    static void addDrug(auto_ptr<LSTMDrugType> drug);
    
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
    LSTMDrugType (const scnXml::PKPDDrug& drugData, uint32_t& bit_start);
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
    
};

} }
#endif
