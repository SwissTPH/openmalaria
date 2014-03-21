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

#ifndef Hmod_WH_Path_Submodels
#define Hmod_WH_Path_Submodels

#include "Global.h"
#include "WithinHost/Pathogenesis/PathogenesisModel.h"

using namespace std;

namespace OM { namespace WithinHost { namespace Pathogenesis {

/// Müller presentation model.
class MuellerPathogenesis : public PathogenesisModel {
public:
  MuellerPathogenesis(double cF) :
    PathogenesisModel(cF) {}
  ~MuellerPathogenesis() {}
  
  virtual double getPEpisode(double timeStepMaxDensity, double totalDensity);

  // Static:
  static void init( const Parameters& parameters );

private:
  ///@brief static vars set by init()
  //@{
  static double rateMultiplier_31;
  static double densityExponent_32;
  //@}
};

/// Pyrogenic threshold presentation model.
class PyrogenPathogenesis : public PathogenesisModel {
protected:
  //!critical density for fever (clinical episodes)
  double _pyrogenThres;
   //! Determine the current pyrogenic threshold.
    virtual void updatePyrogenThres(double totalDensity);

public:
  PyrogenPathogenesis(double cF);
  virtual ~PyrogenPathogenesis() {}
  virtual void summarize (const Host::Human& human);
  virtual double getPEpisode(double timeStepMaxDensity, double totalDensity);
  
  // Static:
  static void init( const OM::Parameters& parameters );
  
protected:
    
    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);
    
private:
  ///@brief Static vars set by init ()
  //@{
  // Ystar2: critical value in determining increase in pyrogenic threshold
  static double Ystar2_13;
  //alpha: factor determining increase in pyrogenic threshold
  static double alpha14;
  //Ystar1: critical value of parasite density in determing increase in pyrog t
  static double Ystar1_26;
  static double smuY;
  //Pyrogenic threshold at birth (Y*0)
  static double initPyroThres;
  //@}
};

/// Predetermined episodes presentation model.
class PredetPathogenesis : public PyrogenPathogenesis {
public:
  PredetPathogenesis (double cF) :
    PyrogenPathogenesis(cF) {}
  ~PredetPathogenesis() {}
  
  virtual double getPEpisode(double timeStepMaxDensity, double totalDensity);
};

} } }
#endif
