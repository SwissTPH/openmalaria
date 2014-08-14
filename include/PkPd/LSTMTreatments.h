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

#include "Global.h"

#include <string>
#include <deque>
#include <vector>
#include <set>
#include <cassert>
#include <memory>
#include <map>
#include <boost/unordered_set.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

using namespace std;

namespace scnXml{
    class Treatments;
    class PKPDMedication;
}
namespace OM { namespace PkPd {
    
/** Information about treatment courses: schedules and dosage tables.
 * 
 * A schedule is a list of times to administrer a drug, along with base drug
 * doses. A dosage table is used to select a multiplier for this base dose
 * according to age (it is assumed for now that age is the only basis on which
 * dosages are selected). */
class LSTMTreatments {
public:
    /** Load treatment data from XML.
     * 
     * Load drug data (LSTMDrugType::init()) first. */
    static void init (const scnXml::Treatments& data);
    
private:
    LSTMTreatments();  // not constructible
};

// for internal use only
struct MedicateData {
    MedicateData () :
        drug(0),
        qty(numeric_limits< double >::signaling_NaN()),
        time(numeric_limits< double >::signaling_NaN()),
        duration(numeric_limits< double >::quiet_NaN())
    {}
    
private:
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        drug & stream;
        qty & stream;
        time & stream;
        duration & stream;
    }
    
    void load( const scnXml::PKPDMedication& med );
    
    size_t drug;      /// Drug type index
    double qty;         /// Quantity of drug prescribed (mg?)
    double time;        /// Time to medicate at (days from start of timestep, may be >= 1 (not this timestep))
    double duration;    /// Duration for IV purposes (use IV admin if a number, oral if is NaN)
    
    friend class Schedule;
    friend class LSTMMedications;
};

/** A class representing PkPd medications allocated but not yet taken (i.e.
 * future doses). Poor adherence is not modeled here, but may be modeled by
 * reducing the number "allocated". */
struct LSTMMedications {
    /// Are there any pending medications?
    inline bool haveMedications(){ return !medicateQueue.empty(); }
    
    /// Call if there are any pending medications.
    void doUpdate( double bodyMass );
    
    //FIXME: call this
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        medicateQueue & stream;
    }
    
private:
    /// All pending medications
    list<MedicateData> medicateQueue;
};

} }
#endif
