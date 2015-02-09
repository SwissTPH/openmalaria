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

#ifndef Hmod_Global
#error "Please include Global.h not SimTime.h directly."
// because checkpoint.h needs this, and we need checkpoint.h
#endif

#include "util/mod.h"
#include "util/checkpoint.h"

#include <iostream>
#include <cassert>
#include <string>

class UnittestUtil;

namespace scnXml {
    class Scenario;
}

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
    
    /// Get length of time in days. Currently this is simple no-op get.
    inline int inDays() const{ return d; }
    
    /// Convert to years
    inline double inYears() const{ return d * (1.0 / DAYS_IN_YEAR); }
    
    /// Convert to time steps (rounding down)
    inline int inSteps() const{ return d / interval; }
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
    
    /** Return this time in time steps modulo some positive integer. */
    inline int moduloSteps(int denominator){
        return util::mod_nn(d / interval, denominator);
    }
    
    /** Return this time in time steps modulo some positive integer. */
    inline int moduloYearSteps(){
        return util::mod_nn(d / interval, steps_per_year);
    }
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        using namespace OM::util::checkpoint;
        d & stream;
    }
    
private:
    int d;      // time in days
    
    // Global constants
    static int interval;        // days per time step
    static size_t steps_per_year;
    static double years_per_step;
    
    friend SimTime mod_nn( const SimTime, const SimTime );
    friend class sim;
};

inline SimTime mod_nn( const SimTime lhs, const SimTime rhs ){
    return SimTime(util::mod_nn(lhs.d, rhs.d));
}
inline int floorToInt( double x ){
	return static_cast<int>(std::floor(x));
}

/** Encapsulation of SimTime static members. */
class sim {
public:
    ///@brief Accessors, all returning a copy to make read-only
    //@{
    /** Time at the beginning of a time step update.
     *
     * This is what is mostly used during an update. It is never negative and
     * increases throughout the simulation. */
    static inline SimTime ts0(){
        assert(in_update);      // should only be used during updates
        return time0;
    }
    /** Time at the end of a time step update.
     * 
     * During an update, ts0() + oneTS() = ts1(). Neither this nor ts0 should
     * be used outside of updates. */
    static inline SimTime ts1(){
        assert(in_update);      // should only be used during updates
        return time1;
    }
    /** Time during init, monitoring and intervention deployment (i.e. whenever
     * a human or mosquito update is not in progress).
     *
     * This is equal to ts1() from the last update and ts0() from this update.
     */
    static inline SimTime now(){
        assert(!in_update);     // only for use outside of step updates
        return time0;   // which is equal to time1 outside of updates, but that's a detail
    }
    /** now() - 1 with same requirement only for use outside of updates. */
    static inline SimTime prevNow(){ return now() - oneTS(); }
    /** During updates, this is ts0; between, this is now. */
    static inline SimTime nowOrTs0(){ return time0; }
    /** During updates, this is ts1; between, this is now. */
    static inline SimTime nowOrTs1(){ return time1; }
    /** During updates, this is ts0; between, it is now - 1. */
    static inline SimTime latestTs0(){ return time1 - oneTS(); }
    
    /** Time relative to the intervention period. Some events are defined
     * relative to this time rather than simulation time, and since the
     * difference is not known until after the warmup period of the simulation
     * it is easier to track the two separately.
     * 
     * This is a large negative number until the intervention period starts,
     * at which time it jumps to zero. */
    static inline SimTime intervNow(){ return interv_time; }
    
    /** One time step (currently either one or five days). */
    static inline SimTime oneTS(){ return SimTime(SimTime::interval); }
    
    /** The number of time steps in one year. */
    static inline size_t stepsPerYear(){ return SimTime::steps_per_year; }
    
    /** A cached value: one year divided by one time step. */
    static inline double yearsPerStep(){ return SimTime::years_per_step; }
    
    static inline SimTime maxHumanAge(){ return max_human_age; }
    //@}
    
    ///@brief Constructors, for convenience
    //@{
    /** Duration zero. */
    static inline SimTime zero(){ return SimTime(0); }
    
    /** One day. */
    static inline SimTime oneDay(){ return SimTime(1); }
    
    /** One year. See SimTime::DAYS_IN_YEAR. */
    static inline SimTime oneYear(){ return SimTime(SimTime::DAYS_IN_YEAR); }
    
    /** Special value representing a time point always in the past, such that
     * never() + x < zero() and x - never() will not to overflow for all valid
     * simulation times x (including any value now() may take as well as
     * never() and future()). */
    static inline SimTime never(){ return SimTime(); }
    
    /** Special value representing a time point always in the future, such that
     * now() < future() and now() + future() does not overflow. */
    static inline SimTime future(){ return SimTime(0x3FFFFFFF); }
    
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
    
    /** Convert from years, rounding down to the next time step. */
    static inline SimTime fromYearsD(double years){
        return fromTS( floorToInt(stepsPerYear() * years) );
    }
    
    /** Convert. */
    static inline SimTime fromTS(int ts){ return oneTS() * ts; }
    
    /** Round to the nearest time-step, where input is in days. */
    static inline SimTime roundToTSFromDays(double days){
        return fromTS( floorToInt( days / SimTime::interval + 0.5 ));
    }
    //@}
    
    ///@brief Conversion functions, for convenience
    //@{
    /** Convert some number of days to some number of time steps (integer
     * division). */
    static inline int daysToSteps(int days){ return days / SimTime::interval; }
    //@}
    
private:
    static void init( const scnXml::Scenario& scenario );
    
    static SimTime max_human_age;   // constant
    // Global variables
#ifndef NDEBUG
    static bool in_update;       // only true during human/population/transmission update
#endif
    static SimTime time0;
    static SimTime time1;
    static SimTime interv_time;
    
    friend class Simulator;
    friend class ::UnittestUtil;
};

}
#endif
