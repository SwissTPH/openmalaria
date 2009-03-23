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

#ifndef Hmod_morbidity
#define Hmod_morbidity

#include <iostream>

using namespace std;

namespace Morbidity {
  /** Types of infection; correspond roughly to those in doCM.
   * 
   * The following are flags:
   * @enum NON_MALARIA = non-malaria infection
   * @enum MALARIA = malaria infection
   * @enum INDIRECT_MORTALITY = Death caused by indirect effects of malaria
   * @enum COMPLICATED = Severe malaria or a coinfection
   * 
   * The following are possible output values:
   * @enum NONE = no infection
   * @enum NON_MALARIA = non-malaria infection
   * 
   * The following output values (malaria infections) may additionally have
   * flag INDIRECT_MORTALITY set:
   * @enum UNCOMPLICATED
   * @enum SEVERE
   * @enum COINFECTION
   */
  enum Infection {
    NONE		= 0,
    
    NON_MALARIA		= 0x1,
    MALARIA		= 0x2,
    
    // Flags indicating morbidity severity:
    INDIRECT_MORTALITY	= 0x4,
    COMPLICATED		= 0x8,
    
    UNCOMPLICATED	= MALARIA | 0x10,
    
    SEVERE		= MALARIA | COMPLICATED | 0x10,
    COINFECTION		= MALARIA | COMPLICATED | 0x20,
  };
}

/*! Morbidity Model abstract base class. */
class MorbidityModel {
public:
  // static:
  /// Calls static init on all MorbidityModels.
  static void initModels();
  
  /** Create a sub-class instance, dependant on global options.
   * 
   * @param cF = Comorbidity factor (currently set in Human). */
  static MorbidityModel* createMorbidityModel(double cF);
  
  // non-static
  MorbidityModel(double cF);
  Morbidity::Infection infectionEvent(double ageYears, double totalDensity, double timeStepMaxDensity);
  bool indirectDeath(double ageYears);
  virtual double getPyrogenThres();
  virtual void write(ostream& out) const;
  virtual void read(istream& in);
  
private:	// static
  //comorbidity prevalence at birth as a risk factor for indirect
  static double indirRiskCoFactor_18;
  //sevMal: critical density for severe malaria episode (Y*B1)
  static double sevMal_21;
  //Critical age for co-morbidity (for both severe and indirect)
  static double critAgeComorb_30;
  //comorbidity prevalence at birth as a risk factor for severe
  static double comorbintercept_24;
  
protected:	// non-static
  virtual double getPEpisode(double timeStepMaxDensity, double totalDensity)=0;
  
  //! comorbidity factor for heterogeneity 
  double _comorbidityFactor; 
};

#endif
