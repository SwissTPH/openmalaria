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

#ifndef Hmod_LSTMPkPdModel
#define Hmod_LSTMPkPdModel

#include "PkPd/Drug/LSTMDrug.h"

#include <boost/ptr_container/ptr_vector.hpp>

namespace scnXml{
    class PKPDMedication;
}
namespace OM {
namespace Host{
    class Human;
}
namespace PkPd {
using boost::ptr_vector;

struct MedicateData {
    MedicateData () :
        drug(0),
        qty(numeric_limits< double >::signaling_NaN()),
        time(numeric_limits< double >::signaling_NaN()),
        duration(numeric_limits< double >::quiet_NaN())
    {}
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        drug & stream;
        qty & stream;
        time & stream;
        duration & stream;
    }
    
private:
    void load( const scnXml::PKPDMedication& med );
    
    inline MedicateData multiplied( double doseMult ){
        MedicateData r( *this );
        r.qty *= doseMult;
        return r;
    }
    
    size_t drug;      /// Drug type index
    double qty;         /// Quantity of drug prescribed (mg when oral, mg/kg when IV)
    double time;        /// Time to medicate at, in days (0 means start of time step, may be >= 1 (thus not today))
    double duration;    /// Duration for IV purposes, in days (use 0 or NaN to indicate oral dose)
    
    friend struct Schedule;
    friend class LSTMModel;
    friend class ::UnittestUtil;
};

/** Pharmacokinetic and pharmacodynamics interface, used by each human's
 * within-host model.
 *
 * This class holds per-human data: prescribed medications and drugs in the
 * body.
 * 
 * Calling order each day:
 *  * prescribe()
 *  * medicate()
 *  * getDrugFactor() for each infection
 *  * decayDrugs()
 */
class LSTMModel {
public:
    /// Static initialisation
    static void init ( const scnXml::Scenario& scenario );
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        checkpoint (stream);
    }
    
    /** Prescribe a patient a course of drugs, via the Pk/Pd model
     * 
     * Note: doses sizes are modified according to age via the dosage
     * table given at the time this function is called.
     *
     * @param schedule Index of a treatment schedule
     * @param dosages Index of a dosage table
     * @param age Age of human at start of time step in years
     * @param body_mass Mass of human in kg
     */
    void prescribe(size_t schedule, size_t dosages, double age, double body_mass);
    
    /** Medicate drugs: human takes prescribed drugs which are to be taken this
     * day.
     * 
     * @param body_mass Mass of human in kg
     * 
     * Note: poor adherence on the part of the patient is not modeled here; to
     * model, prescribe with a "poor adherence" schedule.
     */
    void medicate(double body_mass);
    
    /** Get concentration of the drug at the beginning of the day. */
    double getDrugConc (size_t drug_index) const;
    
    /** This is how drugs act on infections.
     *
     * Each time step, on each infection, the parasite density is multiplied by
     * the return value of this infection. The WithinHostModels are responsible
     * for clearing infections once the parasite density is negligible. */
    double getDrugFactor (uint32_t genotype);
    
    /** After any resident infections have been reduced by getDrugFactor(),
     * this function is called to update drug levels to their effective level
     * at the end of the day, as well as clear data once drug concentrations
     * become negligible. */
    void decayDrugs ();
    
    /** Make summaries of drug concentration data. */
    void summarize( const Host::Human& human ) const;
    
private:
    /** Medicate drugs to an individual, which act on infections the following
     * time steps, until rendered ineffective by decayDrugs().
     *
     * \param typeIndex The index of drug type data (what LSTMDrugType::findDrug() returns).
     * \param qty The quantity in either mg (if oral dose) or mg/kg (if IV).
     * \param time Time in days since start of this time step to medicate at
     * \param duration  Duration in days. 0 or an NaN indicates no duration.
     * \param bodyMass Weight of human in kg
     * 
     * Due to the fact we're using a discrete time step model, the case-management
     * update (calling medicate) and within-host model update (calling
     * getDrugFactor) cannot [easily] have immediate effects on each other. The
     * implementation we use is that the within-host model update (calculating
     * new infection densities) happens first; hence medicate() will always be
     * called after getDrugFactor in a time step, and a time of zero means the
     * dose has effect from the start of the following time step. */
    void medicateDrug(size_t typeIndex, double qty, double time, double duration, double bodyMass);
  
    void checkpoint (istream& stream);
    void checkpoint (ostream& stream);
    
    typedef ptr_vector<LSTMDrug> DrugVec;
    /// Drugs with non-zero blood concentrations:
    DrugVec m_drugs;
    
    /// All pending medications
    list<MedicateData> medicateQueue;
    
    friend class ::UnittestUtil;
};

} }
#endif