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

#ifndef Hmod_LSTMDrug
#define Hmod_LSTMDrug

#include "Global.h"

using namespace std;

namespace OM {

namespace PkPd {
struct DoseParams;
}
namespace util {
namespace checkpoint {
void operator& (multimap<double,PkPd::DoseParams> x, ostream& stream);
void operator& (multimap<double,PkPd::DoseParams>& x, istream& stream);
}
}

namespace PkPd {

/** Describes an oral or IV dose.
 *
 * If duration > 0, describes an IV dose. qty refers to the infusion rate
 * in mg/kg/day.
 *
 * If duration = 0, describes an oral dose. qty refers to the concentration
 * in mg/l.
 */
struct DoseParams {
    DoseParams() : qty(0.0), duration(0.0) {}
    DoseParams( double r, double d ): qty(r), duration(d) {}
    double qty;               // infusion rate (mg/kg/day) or dose size (mg/l)
    //NOTE: duration is now only needed when checking for overlapping doses
    // and not within factor/concentration calculation code.
    double duration;  // units: days

    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        qty & stream;
        duration & stream;
    }
};

/** A class holding pkpd drug use info.
 *
 * Each human has an instance for each type of drug present in their blood. */
class LSTMDrug {
public:
    /** Create a new instance. */
    LSTMDrug ();
    
    /** Indicate a new medication this time step.
     *
     * Converts qty in mg to concentration, and stores along with time (delay past
     * the start of the current time step) in the doses container.
     * 
     * @param time Time of administration, in days (should be at least 0 and
     *  less than 1).
     * @param qty Amount of active ingredient, in mg
     * @param volDist Volume of distribution in l
     */
    void _medicate (double time, double qty, double volDist);
    /** Indicate a new medication via IV this time step.
     *
     * @param time Time of start of administration, in days (should be at least
     *  0 and less than 1 (although time+duration may be greater than 1).
     * @param duration Duration of IV, in days. Drug is assumed to be
     *  administered at a constant rate of this duration.
     * @param qty Quantity of active ingredient, in mg/kg
     */
    void medicateIV (double time, double duration, double qty);

    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        doses & stream;
        checkpoint(stream);
    }

protected:
    virtual void checkpoint (istream& stream){}
    virtual void checkpoint (ostream& stream){}
    
    typedef multimap<double,DoseParams> DoseMap;

    /** Check whether an IV dose needs to be split into multiple doses
     * (over two days or when an oral dose occurs in the middle).
     * If necessary, split.
     *
     * Only check against lastInserted. */
    void check_split_IV( DoseMap::iterator lastInserted );


    /** List of each dose given today (and possibly tomorrow), ordered by time.
     * First parameter (key) is time in days, second describes dose.
     *
     * Used in calculateDrugFactor temporarily,
     * and in updateConcentration() to update concentration. */
    DoseMap doses;
};

}
}
#endif
