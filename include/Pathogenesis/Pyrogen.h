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

#ifndef Hmod_pyrogenMorb
#define Hmod_pyrogenMorb

#include "Global.h"
#include "Pathogenesis/PathogenesisModel.h"

using namespace std;

namespace OM { namespace Pathogenesis {

/*! Pyrogenic threshold presentation model.
*/
class PyrogenPathogenesis : public PathogenesisModel {
protected:
  //!critical density for fever (clinical episodes)
  double _pyrogenThres;
   //! Determine the current pyrogenic threshold.
    virtual void updatePyrogenThres(double totalDensity);

public:
  PyrogenPathogenesis(double cF);
  virtual ~PyrogenPathogenesis() {}
  virtual void summarize (Survey& survey, SurveyAgeGroup ageGroup);
  virtual double getPEpisode(double timeStepMaxDensity, double totalDensity);
  
  // Static:
  static void init();
  
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

} }
#endif