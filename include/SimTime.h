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

class UnittestUtil;

namespace OM {

/******************************************************************************
 * Class encapsulating simulation time (as in days and dates, not time-of-day).
 * 
 * Time steps, days and dates are derived from this. The values and units of
 * internal variables are an implementation detail (i.e. code outside this
 * class should not need to know).
 * 
 * Type represents relative times (durations) and absolute times (duration
 * since start of the simulation or since the start of the intervention
 * period).
 *****************************************************************************/
class SimTime {
    /** Construct, from a time in days. */
    explicit SimTime( int days ) : d(days) {}
    
public:
    /// Number of days in a year; defined as 365 (leap years are not simulated).
    enum { DAYS_IN_YEAR = 365 };
    
    /** Default construction; same as sim::never(). */
    SimTime() : d(-0x3FFFFFFF) {}
    
    ///@brief Conversions to other types/units
    //NOTE: these methods provide good documentation of the types of things
    //one does with SimTimes (besides comparing with other SimTimes).
    //@{
    /// Get raw value (currently days; not guaranteed not to change). Same value as checkpointed.
    inline int raw() const{ return d; }
    
    /// Convert to TimeStep
    inline util::TimeStep ts() const{ return util::TimeStep::fromDays(d); }
    
    /// Get length of time in days. Currently this is simple no-op get.
    inline int inDays() const{ return d; }
    
    /// Convert to years
    inline double inYears() const{ return d * (1.0 / DAYS_IN_YEAR); }
    
    /// Get array index in time steps (for dense arrays, involves conversion)
    inline int32_t indexTS() const{
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
    // scale by an integer
    inline SimTime operator*( int scalar )const {
        return SimTime( d * scalar );
    }
    // scale by a double, rounding to nearest
    inline SimTime operator*( double scalar )const {
        return SimTime( static_cast<int>(d * scalar + 0.5) );
    }
    // Divide by another SimTime; result is unitless. Note integer division.
    inline int operator/( const SimTime rhs )const{
        return d / rhs.d;
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
    int d;      // time in days
    
    friend std::ostream& operator<<( std::ostream&, const SimTime );
    friend SimTime mod_nn( const SimTime, const SimTime );
    friend class sim;
};

inline std::ostream& operator<<( std::ostream& stream, const SimTime time ){
    return( stream << time.d );
}
inline SimTime mod_nn( const SimTime lhs, const SimTime rhs ){
    return SimTime(util::mod_nn(lhs.d, rhs.d));
}

/** Encapsulation of SimTime static members. */
class sim {
public:
    ///@brief Accessors, all returning a copy to make read-only
    //@{
    /** Get the time now (i.e. duration since start of simulation, including
     * initialisation period). The following is always true: now() >= zero().
     */
    static inline SimTime now(){ return sim_time; }
    
    static inline SimTime maxHumanAge(){
        return fromTS(util::TimeStep::maxAgeIntervals);
    }
    //@}
    
    ///@brief Constructors, for convenience
    //@{
    /** Duration zero. */
    static inline SimTime zero(){ return SimTime(0); }
    
    /** One day. */
    static inline SimTime oneDay(){ return SimTime(1); }
    
    /** One time step (currently either one or five days). */
    static inline SimTime oneTS(){ return one_step; }
    
    /** One year. See SimTime::DAYS_IN_YEAR. */
    static inline SimTime oneYear(){ return SimTime(SimTime::DAYS_IN_YEAR); }
    
    /** Special value representing a time point always in the past, such that
     * never() + x < zero() and x - never() will not to overflow for all valid
     * simulation times x (including any value now() may take as well as
     * never() and future()). */
    static inline SimTime never(){ return SimTime(); }
    
    /** Duration in days. Should be fast (currently no conversion required). */
    static inline SimTime fromDays(int days){ return SimTime(days); }
    
    /** Convert from a whole number of years. */
    static inline SimTime fromYearsI(int years){
        return SimTime(SimTime::DAYS_IN_YEAR * years);
    }
    
    /** Convert from years to nearest time step. */
    static inline SimTime fromYearsN(double years){
        return roundToTSFromDays(SimTime::DAYS_IN_YEAR * years);
    }
    
    /** Convert. */
    static inline SimTime fromTS(const util::TimeStep ts){
        return SimTime(ts.inDays());
    }
    
    /** Convert. */
    static inline SimTime fromTS(int ts){ return oneTS() * ts; }
    
    /** Round to the nearest time-step, where input is in days. */
    static inline SimTime roundToTSFromDays(double days){
        return fromTS(std::floor( days / one_step.d + 0.5 ));
    }
    //@}
    
private:
    static SimTime sim_time;
    static SimTime one_step;
    
    friend class Simulator;
    friend class ::UnittestUtil;
};

}
#endif
