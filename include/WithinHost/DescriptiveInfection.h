/*
 This file is part of OpenMalaria.

 Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine

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

#ifndef Hmod_DescriptiveInfection
#define Hmod_DescriptiveInfection
#include <fcntl.h>
#include "intervention.h"
#include "WithinHost/Infection.h"
#include "Simulation.h"
#include <fstream>


//Max duration of an infection in intervals. TODO: Consequences for non-5day interval simulations?
const int maxDur=84;

//The maximum parasite density we allow per DescriptiveInfection. Higher values are set to maxDens.
const double maxDens=2000000;


//!  Models of infection.
/*!
  Models related to the within-host dynamics of infections.
*/
class DescriptiveInfection : public Infection {
//TODO: should be private, and immune decay and immune proxies need to be discussed in light of new WIH-models 
public:
  ///@name Static init/cleanup
  //@{
  //! Init constants common to all Phase A (AJTMH 75(2)) infections.
  /*!
    Init constants common to all infections modeled via the original (AJTMH 75(2)
    empirical model.  Using this model, the time step remains 5 days.  Where the simulation
    time step is shorter than 5 days the parasite densities are looked up by rounding down
    to the previous 5 days.

    Once constants are initialised then cumulative distributions of parasite densities and durations of patency
    from the malariatherapy data and also the category boundaries for the grouping
    of time since first positive slide.
  */
  static void initParameters();
  static void clearParameters();
  //@}
  
  ///@name CTOR & DTOR
  //@{
  //! Constructor
  DescriptiveInfection(int simulationTime);
  
  /** Checkpoint-reading constructor */
  DescriptiveInfection (istream& in);
  
  /** Destructor */
  virtual ~DescriptiveInfection();
  //@}
  
  virtual void write (ostream& out) const;
  
  //! Get the last timestep before the infection is cleared.
  /*!
    \return The interval before clearance.
  */
  bool expired () {
    return Simulation::simulationTime >= _startdate+_duration;
  }
  
  //! Get the density of the infection
  double getDensity() { return _density; }

  //! Start date of the infection
  int getStartDate() { return _startdate; }
  
  /** Determines parasite density of an individual infection (5-day timestep update)
   *
   * @param ageInYears Age of human
   * @param cumulativeh Cumulative number of infections
   * @param cumulativeY Previous exposure (cumulative parasite density)
   * @param timeStepMaxDensity (In-out param) Used to return the maximum parasite density over a 5-day interval.
   * @param expInnateImm Density multiplier for innate immunity.
   * @param BSVEfficacy Density multiplier for Blood-Stage Vaccine effect.
   */
  void determineDensities(double ageInYears, double cumulativeh, double cumulativeY, double &timeStepMaxDensity, double expInnateImm, double BSVEfficacy);
  /// Final part of determineDensities calculation.
  /// Separated out to preserve order of operations.
  void determineDensityFinal ();
  
  //! Initialises infection duration.
  /*! 
    Initialises infection duration sampling from log normal distribution using parameters for 53 patients from Georgia.
    Mean log duration of an infection values from AJTM p.9 eq.5.
    \return The duration in simulation intervals.
  */
  int infectionDuration();

  
  //! Write an infection to the checkpointing file.
  /*!
    \param funit Checkpoint file.
  */
  void writeInfectionToFile (fstream& funit);
  
  /// Includes the effect of attenuated infections by SP concentrations, when using IPT
  virtual void IPTattenuateAsexualDensity () {}
  
  /// pre-erythrocytic latent period, in time steps
  //Note: here for convenience; used by DescriptiveInfection and DescriptiveIPT
  static int latentp;
  
protected:
  //! Arbitrary maximum duration of the infection, in timesteps
  int _duration; 
  
private:
  // -----  static  -----
  
  //! Density distributions
  /*!
    Mean Log Parasite Count at time step i for an infection that lasts j days.
    only about one half of the matrix is initialized. (right upper triangle)
  */
  static double meanLogParasiteCount[maxDur*maxDur];
 
  static double sigma0sq; //!< Sigma0^2 in AJTM p.9 eq. 13

  static double xNuStar; //!< XNuStar in AJTM p.9 eq. 13
};

#endif

