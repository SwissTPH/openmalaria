/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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

/** The Global module is only for items wanted nearly everywhere.
 * It is not for miscelaneous stuff and should only include headers wanted nearly everywhere.
 */
#ifndef Hmod_Global
#define Hmod_Global

#ifdef _MSC_VER
#pragma warning(disable: 4290)	// disable some warnings on MSVC
#define finite(x) _finite(x)
#endif

#include "util/checkpoint.hpp"
using namespace OM::util::checkpoint;

namespace OM {
    /// Global stuff. (Would have been a class but can't initialize constants in a class; better
    /// that these values are compile-time known.)
    namespace Global {
	/** Sets parameters in Global and performs some checks. */
	void init ();
	
	///@brief Compile-time constants
	//@{
	/** Value used as the timestep for an event which has never happened.
	*
	* For any simulation timestep, we must have:
	* ( TIMESTEP_NEVER + simulationTime < 0 )
	* but since (x - TIMESTEP_NEVER >= y) is often checked, x - TIMESTEP_NEVER
	* must not overflow for any timestep x (int represents down to -0x7FFFFFFF).
	*/
	const int TIMESTEP_NEVER = -0x3FFFFFFF;
	
	/// Days in a year. Should be a multiple of interval.
	const int DAYS_IN_YEAR = 365;
	//@}
	
	///@brief Variables read from XML configuration which remain constant after set up
	//@{
	/// temporal resolution of simulation, in days
	extern int interval;
	/// Number of timesteps in 5 days (i.e. 5/interval)
	extern int intervalsPer5Days;
	/// Simulation time steps per year (DAYS_IN_YEAR / interval)
	extern size_t intervalsPerYear;
	/// 1.0/intervalsPerYear
	extern double yearsPerInterval;
	/// Maximum age of individuals in a scenario in time intervals
	extern int maxAgeIntervals;
	//@}
	
	/** @brief Simulation variables, made global to save continuously passing.
	 *
	 * Don't add many variables like this.
	 * 
	 * These should be considered read-only outside Simulation. Checkpointed. */
	//@{
	/// Simulation timestep, starting from beginning of initialization.
	/// (Each step is "interval" long).
	extern int simulationTime;
	
	/** Timestep counter during the main simulation.
	*
	* This is <0 during initialization and is incremented from 0 from the start of the intervention
	* period. */
	extern int timeStep;
	//@}
    };
}
#endif
