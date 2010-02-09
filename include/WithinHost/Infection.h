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

#ifndef Hmod_Infection
#define Hmod_Infection

#include "Global.h"
#include "PkPd/Proteome.h"

namespace OM { namespace WithinHost {
    
class Infection {
public:
  static void init();
  
  Infection (uint32_t protID) :
    proteome_ID(protID),
    _density(0.0),
    _cumulativeExposureJ(0.0)
    {}
  Infection (istream& stream);
  virtual ~Infection () {}
  
  //! Get proteome
  inline uint32_t get_proteome_ID() const {
    return proteome_ID;
  }
  
  
  /** @returns A multiplier describing the proportion of parasites surviving
   * immunity effects this timestep.
   * 
   * Note that in the Descriptive model this multiplies log(density), but the
   * new density has no effect on future densities, wheras the Empirical model
   * multiplies the actual density (which then affects density on the following
   * timestep). */
  double immunitySurvivalFactor (double ageInYears, double cumulativeh, double cumulativeY);
  
  /// Checkpointing
  template<class S>
  void operator& (S& stream) {
      checkpoint (stream);
  }
  
protected:
    virtual void checkpoint (ostream& stream);
    
  //! Proteome (used in a different situation than genotype) 
  uint32_t proteome_ID;
  
  //! Current density of the infection
  double _density;
  
  //! Cumulative parasite density, since start of this infection
  double _cumulativeExposureJ;
  
  /// @brief Static data set by init
  //@{
public:
  static float cumulativeYstar; //!< Critical value for immunity trigger (cumulative densities)
  static float cumulativeHstar; //!< Critical value for immunity trigger (cumulative inoculations)
  
private:
  static double alpha_m; //!< Maternal protection at birth
  
  /*!
  More or less (up to 0.693) inverse quantity of alphaMStar (AJTM p. 9 eq. 12),
  decay rate of maternal protection in years^(-1).
  */
  static double decayM;
  //@}
};

} }
#endif
