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
#include "Infection.h"
#include <fstream>

typedef struct genotype genotype;

struct genotype {
   int ID;
//!In order to save memory, we just define the ID of the genotype. Attributes of the
//!genotype can be accessed via arrays in mod_intervention.
//!(e.g. freq = mod_intervention.GenotypeFreq(iTemp%iData%gType%ID)
//!attributes are:
//!freq: Probability of being infected by this specific genotype
//!ACR: Probability of being cured (due to SP)
//!proph: Prophylactic effect of SP (measured in time steps)
//!tolperiod: time window of tolerance period
//!SPattenuation: Factor of how parasites are attenuated  by SP (genotype specific)
};


//!  Models of infection.
/*!
  Models related to the within-host dynamics of infections.
*/
class DescriptiveInfection : public Infection {
  //! Proteome (used in a different situation than genotype) 
  ProteomeInstance* _proteome; 

 //TODO: should be private, and immune decay and immune proxies need to be discussed in light of new WIH-models 
  public:

   DescriptiveInfection();
  ~DescriptiveInfection();

  //! Constructor
  /*!
    \param lastSPdose Time interval of last SP Dose.
  */
  DescriptiveInfection(int lastSPdose, int simulationTime);
  
  void write (ostream& out);
  void read (istream& in);
  
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
  
  //! Get the last timestep before the infection is cleared.
  /*!
    \return The interval before clearance.
  */
  int getEndDate();

  //! Get the density of the infection
  double getDensity() { return _density; };

  //! Get proteome
  ProteomeInstance* getProteome() const;

  bool getSPattenuate() { return _SPattenuate; };
  
  void setDensity(double density) { _density = density;};

  int getGenoTypeID() { return _gType.ID; };

  //! Start date of the infection
  int getStartDate() { return _startdate; };
  
  int getDuration() { return _duration; };

  double getCumulativeExposureJ() {return _cumulativeExposureJ;};

  void setCumulativeExposureJ(double exposure) { _cumulativeExposureJ=exposure; };

  //! Determines parasite density of an individual infection.
  /*!
    \param cumulativeY Previous exposure, in cumulative number of parasites.
    \param ageyrs Age in years.
    \param cumulativeH cumulative number of inoculations (inoculation equals infection ?)

  */
  void determineDensities(int simulationTime, double cumulativeY, double ageyears, double cumulativeh, double *timeStepMaxDensity);

  /*
  //! Dummy function. Should override \sa determineDensities.
  */ 
  double determineWithinHostDensity();

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

  double getAlpha_m() const {return alpha_m;}
  double getDecayM() const {return decayM;}
  double getSigma0sq() const {return sigma0sq;}
  double getXNuStar() const {return xNuStar;}
  double getMeanLogParasiteCount(int pos) const {return meanLogParasiteCount[pos];}
  float getCumulativeHstar() const {return cumulativeHstar;};
  float getCumulativeYstar() const {return cumulativeYstar;};


 private:

  //! Density distributions
  /*!
    Mean Log Parasite Count at time step i for an infection that lasts j days.
    only about one half of the matrix is initialized. (right upper triangle)
  */
  static double meanLogParasiteCount[maxDur*maxDur];
 
  static double alpha_m; //!< Maternal protection at birth

  /*!
    More or less (up to 0.693) inverse quantity of alphaMStar (AJTM p. 9 eq. 12),
    decay rate of maternal protection in years^(-1).
  */
  static double decayM;

  static double sigma0sq; //!< Sigma0^2 in AJTM p.9 eq. 13

  static double xNuStar; //!< XNuStar in AJTM p.9 eq. 13
 
  //! Sampled duration of the infection
  int _duration; 
  //! Start date of the infection
  int _startdate; 
  //! Current density of the infection
  double _density;
  //! Cumulative parasite density, since start of this infection
  double _cumulativeExposureJ; 
  //! Genotype responsible for infection
  genotype _gType;
  //! IPTi parameter (indicator for attenuation).
  bool _SPattenuate; 


};

#endif

