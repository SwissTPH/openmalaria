/*
 This file is part of OpenMalaria.

 Copyright (C) 2005-2012 Swiss Tropical Institute and Liverpool School Of Tropical Medicine

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

#ifndef Hmod_TimeStep
#define Hmod_TimeStep

#include "util/mod.h"

#include <iostream>
#include <cassert>

namespace OM {
namespace util {

/** Class encapsulating a time step to add type safety. Conversion to this
    * class must be explicit; also potentially conversion from this class?
    * 
    * Type represents both times and durations.
    */
class TimeStep {
    int _ts;        // time step
    
    // encapsulates a value-type such that only parent class can set it
    template<typename T>
    class ReadOnly {
        T value;
        // only friends can set value:
        inline void operator=(T n){
            value = n;
        }
        friend class TimeStep;
    public:
        // anyone can get value
        // (uses implicit cast, so no extra syntax needed to get value)
        inline operator T () const{
            return value;
        }
    };
    
public:
    /** Sets parameters and performs some checks.
     * 
     * @param daysPerTimeStep Length of an interval
     * @param maxAgeYears Maximum age of a human
     */
    static void init (int daysPerTimeStep, double maxAgeYears);

    /// Days in a year. Should be a multiple of interval.
    enum { DAYS_IN_YEAR = 365 };
    
    /** Simulation time step.
     * 
     * Set-up of populations occurs at time 0.
     * 
     * Each update occurs between time t-1 and t (where t is
     * TimeStep::simulation), so the value of this variable during the first
     * update is 1. */
    static TimeStep simulation;

    /** Timestep counter during the intervention period of the simulation.
    *
    * This is negative during initialization and is incremented from 0
    * from the start of the intervention period (once this period starts,
    * TimeStep::interventionPeriod is always TimeStep::simulation - x, where
    * x is the time at which the intervention period starts. Note that the value
    * of x is not always known definitively before the intervention period is
    * started.
    * 
    * Surveys happen according to times given in the XML, which is this
    * variable's measure. The same is true for intervention times given in the
    * XML. */
    static TimeStep interventionPeriod;

    ///@brief Variables read from XML configuration which remain constant after set up
    //@{
    /// temporal resolution of simulation, in days
    static ReadOnly<int> interval;
    /// 1.0/intervalsPerYear
    static ReadOnly<double> yearsPerInterval;
    /// Number of timesteps in 5 days (i.e. 5/interval)
    static TimeStep intervalsPer5Days;
    /// Simulation time steps per year (DAYS_IN_YEAR / interval)
    static TimeStep intervalsPerYear;
    /// Maximum age of individuals in a scenario in time intervals
    static TimeStep maxAgeIntervals;
    /// Same as intervalsPerYear, but as an integer: useful for array indecies and lengths
    static ReadOnly<int> stepsPerYear;
    //@}

    /** Special values used as the time for an event which has never happened
     * and an event which will always be in the future.
     *
     * For any simulation timestep, we must have: ( never + simulation < 0 ),
     * but since (x - never >= y) is often checked, x - never must not overflow
     * for any timestep x (int represents down to -0x7FFFFFFF).
     * 
     * We must also always have ( simulation < future ). */
    static const TimeStep never, future;
    
    /** Initialize to never. */
    TimeStep();
    /** Convert an integer to a TimeStep. */
    explicit TimeStep( int ts ) : _ts(ts) {}
    
    /// Convert a real number of timesteps to the nearest TimeStep
    static inline TimeStep fromNearest( double d ){
        return TimeStep( static_cast<int>(d+0.5) );
    }
    /// Convert a number of days to TimeStep type, rounding to nearest.
    static TimeStep fromDaysNearest( double d );
    /// Convert a number of days to TimeStep type, rounding down.
    static TimeStep fromDays( double d );
    /// Convert a whole number of years to TimeStep type
    static TimeStep fromYears( int y ){
        assert( intervalsPerYear._ts != 0 );    // not initialized yet
        return TimeStep( y * intervalsPerYear._ts );
    }
    /** Convert a floating-point number of years to TimeStep type (rounding
     * down where rounding is necessary). */
    static TimeStep fromYears( double y ){
        assert( intervalsPerYear._ts != 0 );    // not initialized yet
        return TimeStep( static_cast<int>(y * intervalsPerYear._ts) );
    }
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        using namespace OM::util::checkpoint;
        _ts & stream;
    }
    
    // self-modifying operators
    inline void operator++() {
        ++_ts;
    }
    inline void operator--() {
        --_ts;
    }
    inline void operator+=( const TimeStep rhs ) {
        _ts += rhs._ts;
    }
    
    // arithmatic modifiers returning a copy
    inline TimeStep operator-()const {
        return TimeStep( -_ts );
    }
    inline TimeStep operator-( const TimeStep rhs )const {
        return TimeStep( _ts - rhs._ts );
    }
    inline TimeStep operator+( const TimeStep rhs )const {
        return TimeStep( _ts + rhs._ts );
    }
    // scale the TimeStep by a double, rounding to nearest
    inline TimeStep operator*( double rhs )const {
        return TimeStep( static_cast<int>(_ts * rhs + 0.5) );
    }
    
    // boolean operators
    inline bool operator==( const TimeStep rhs )const {
        return  _ts == rhs._ts;
    }
    inline bool operator!=( const TimeStep rhs )const {
        return  _ts != rhs._ts;
    }
    inline bool operator>( const TimeStep rhs )const {
        return  _ts > rhs._ts;
    }
    inline bool operator>=( const TimeStep rhs )const {
        return  _ts >= rhs._ts;
    }
    inline bool operator<( const TimeStep rhs )const {
        return  _ts < rhs._ts;
    }
    inline bool operator<=( const TimeStep rhs )const {
        return  _ts <= rhs._ts;
    }
    
    /// Get value in days
    inline int inDays() {
        return _ts * interval;
    }
    /// Get value in years
    inline double inYears() {
        return _ts * yearsPerInterval;
    }
    /// Get value in time steps as an integer (marginally faster than inDays() and inYears())
    inline const int asInt()const {
        return _ts;
    }
    
    friend std::ostream& operator<<( std::ostream&, const TimeStep );
    friend int mod_nn( const TimeStep, int );
    friend int mod( const TimeStep, int );
    friend TimeStep mod_nn( const TimeStep, const TimeStep );
    friend TimeStep mod( const TimeStep, const TimeStep );
};

inline std::ostream& operator<<( std::ostream& stream, const TimeStep ts ){
    return( stream << ts._ts );
}

// returns an int since common usage is to get an array index
inline int mod_nn( const TimeStep ts, int rhs ) {
    return util::mod_nn(ts._ts, rhs);
}
inline int mod( const TimeStep ts, int rhs ) {
    return util::mod(ts._ts, rhs);
}
inline TimeStep mod_nn( const TimeStep lhs, const TimeStep rhs ){
    return TimeStep(util::mod_nn(lhs._ts, rhs._ts));
}
inline TimeStep mod( const TimeStep lhs, const TimeStep rhs ){
    return TimeStep(util::mod(lhs._ts, rhs._ts));
}

}
}
#endif
