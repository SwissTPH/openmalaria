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

#include "global.h"
#include "proteome.h"


//Max duration of an infection in intervals. TODO: Consequences for non-5day interval simulations?
const int maxDur=84;

//The maximum parasite density we allow per DescriptiveInfection. Higher values are set to maxDens.
const double maxDens=2000000;


//TODO
double sampleFromLogNormal (double normp, double meanlog, double stdlog);

class Infection {
public:
  // Used in drug.cpp
  //@{
  //! Get proteome
  virtual ProteomeInstance* getProteome() const =0;
  //@}
  
  // Used in oldWithinHostModel.cpp
  //@{
  //! Start date of the infection
  virtual int getStartDate() =0;
  virtual int getDuration() =0;
  
  virtual double getDensity() =0;
  virtual void setDensity(double density) =0;
  
  virtual double getCumulativeExposureJ() =0;
  virtual void setCumulativeExposureJ(double exposure) =0;
  
  virtual bool getSPattenuate() =0;
  
  virtual int getGenoTypeID() =0;
  
  virtual double determineWithinHostDensity() =0;
  
  virtual double getAlpha_m() const =0;
  virtual double getDecayM() const =0;
  
  virtual double getSigma0sq() const =0;
  virtual double getXNuStar() const =0;
  virtual double getMeanLogParasiteCount(int pos) const =0;
  virtual float getCumulativeHstar() const =0;
  virtual float getCumulativeYstar() const =0;
  //@}
  
  // Used in human.cpp
  //@{
  //! Get the last timestep before the infection is cleared.
  /*!
  \return The interval before clearance.
   */
  virtual int getEndDate() =0;
  
  static float cumulativeYstar; //!< Critical value for immunity trigger (cumulative densities)
  static float cumulativeHstar; //!< Critical value for immunity trigger (cumulative inoculations)
  
  virtual void write (ostream& out) =0;
  virtual  void read (istream& in) =0;
  //@}
};

#endif
