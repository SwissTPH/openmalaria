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

#ifndef Hmod_OM_SimTime
#define Hmod_OM_SimTime

#include "util/TimeStep.h"      // for conversion

#include <iostream>
#include <cassert>

namespace OM {

/******************************************************************************
 * Class encapsulating simulation time.
 * 
 * Time steps, days and dates are derived from this. The values and units of
 * internal variables are an implementation detail (i.e. code outside this
 * class should not need to know).
 * 
 * Type represents both times (from some epoch) and durations.
 *****************************************************************************/
class SimTime {
    /** Construct, from a time in days. */
    explicit SimTime( int days ) : d(days) {}
    
public:
    ///@brief Simple arithmatic modifiers (all return a copy)
    //@{
    inline SimTime operator-()const {
        return SimTime( -d );
    }
    inline SimTime operator-( const SimTime rhs )const {
        return SimTime( d - rhs.d );
    }
    inline SimTime operator+( const SimTime rhs )const {
        return SimTime( d + rhs.d );
    }
    //@}
    
    ///@brief Comparators between two SimTimes (all return a boolean)
    //@{
    inline bool operator==( const SimTime rhs )const {
        return  d == rhs.d;
    }
    inline bool operator!=( const SimTime rhs )const {
        return  d != rhs.d;
    }
    inline bool operator>( const SimTime rhs )const {
        return  d > rhs.d;
    }
    inline bool operator>=( const SimTime rhs )const {
        return  d >= rhs.d;
    }
    inline bool operator<( const SimTime rhs )const {
        return  d < rhs.d;
    }
    inline bool operator<=( const SimTime rhs )const {
        return  d <= rhs.d;
    }
    //@}
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        using namespace OM::util::checkpoint;
        d & stream;
    }
    
private:
    uint32_t d;      // time in days
    
    friend class sim;
};

/** Encapsulation of SimTime static members. */
class sim {
public:
    ///@brief Accessors, all returning a copy to make read-only
    //@{
    /** Get the time now (i.e. duration since start of simulation, including
     * initialisation period). The following is always true: now() >= zero().
     */
    static inline SimTime now(){ return fromTS(util::TimeStep::simulation); }
    //@}
    
    ///@brief Constructors, for convenience
    //@{
    /** Duration zero. */
    static inline SimTime zero(){ return SimTime(0); }
    
    /** Special value representing a time point always in the past, such that
     * the following is always true: never() + now() < zero() . Additionally,
     * x - never() is guaranteed not to overflow for any SimTime x >= 0. */
    static inline SimTime never(){
        //TODO: this implies a maximum value for now()
        return SimTime(-0x3FFFFFFF);
    }
    
    /** Duration in days. Should be fast (currently no conversion required). */
    static inline SimTime fromDays(int days){ return SimTime(days); }
    
    /** Convert. */
    static inline SimTime fromTS(const util::TimeStep ts){ return SimTime(ts.inDays()); }
    //@}
    
private:
//     static SimTime m_now;
};

}
#endif
