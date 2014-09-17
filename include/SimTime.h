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
    explicit SimTime( int32_t days ) : d(days) {}
    
public:
    /** Default construction; same as sim::never(). */
    SimTime() : d(-0x3FFFFFFF) {}
    
    ///@brief Conversions to other types/units
    //@{
    /// Get raw value (currently days; not guaranteed not to change). Same value as checkpointed.
    inline int32_t raw() const{ return d; }
    
    /// Convert to TimeStep
    inline util::TimeStep ts() const{ return util::TimeStep::fromDays(d); }
    
    /// Convert to years
    inline double inYears() const{ return d * (1.0 / 365); }
    
    /// Get array index in time steps (for dense arrays, involves conversion)
    inline size_t indexTS() const{
        return d / util::TimeStep::interval;
    }
    //@}
    
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
    
    ///@brief Self-modifying arithmatic
    //@{
    inline void operator+=( const SimTime rhs ) {
        d += rhs.d;
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
    int32_t d;      // time in days
    
    friend std::ostream& operator<<( std::ostream&, const SimTime );
    friend class sim;
};

inline std::ostream& operator<<( std::ostream& stream, const SimTime time ){
    return( stream << time.d );
}

/** Encapsulation of SimTime static members. */
class sim {
public:
    ///@brief Accessors, all returning a copy to make read-only
    //@{
    /** Get the time now (i.e. duration since start of simulation, including
     * initialisation period). The following is always true: now() >= zero().
     */
    static inline SimTime now(){ return fromTS(util::TimeStep::simulation); }
    
    static inline SimTime maxHumanAge(){ return fromTS(util::TimeStep::maxAgeIntervals); }
    //@}
    
    ///@brief Constructors, for convenience
    //@{
    /** Duration zero. */
    static inline SimTime zero(){ return SimTime(0); }
    
    /** One time step (currently either one or five days). */
    static inline SimTime oneTS(){ return SimTime(util::TimeStep::interval); }
    
    /** One year. This is defined as 365 days in the simulator. */
    static inline SimTime oneYear(){ return SimTime(365); }
    
    /** Special value representing a time point always in the past, such that
     * never() + x < zero() and x - never() will not to overflow for all valid
     * simulation times x (including any value now() may take as well as
     * never() and future()). */
    static inline SimTime never(){
        return SimTime();
    }
    
    /** Duration in days. Should be fast (currently no conversion required). */
    static inline SimTime fromDays(int32_t days){ return SimTime(days); }
    
    /** Convert. */
    static inline SimTime fromTS(const util::TimeStep ts){ return SimTime(ts.inDays()); }
    //@}
    
private:
//     static SimTime m_now;
};

}
#endif
