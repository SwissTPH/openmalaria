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

#ifndef Hmod_PathogenesisModel
#define Hmod_PathogenesisModel

#include <iostream>
#include <vector>

using namespace std;

namespace Pathogenesis {
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
    
    // Flags indicating infection severity:
    INDIRECT_MORTALITY	= 0x4,
    COMPLICATED		= 0x8,
    
    UNCOMPLICATED	= MALARIA | 0x10,
    
    SEVERE		= MALARIA | COMPLICATED | 0x10,
    COINFECTION		= MALARIA | COMPLICATED | 0x20,
  };
}

/*! PathogenesisModel abstract base class. 
 *
 * Previously named MorbidityModel and PresentationModel. */
class PathogenesisModel {
public:
  // static:
  /// Calls static init on correct PathogenesisModel.
  static void init();
  
  /** Create a sub-class instance, dependant on global options.
  * 
  * @param cF = Comorbidity factor (currently set in Human). */
  static PathogenesisModel* createPathogenesisModel(double cF);
  /** Create a sub-class instance, loading from a checkpoint. */
  static PathogenesisModel* createPathogenesisModel(istream& in);

  /** Called for each birth; returns true if infant dies due to mother's
   * infection. */
  static bool eventNeonatalMortality();
  /** Calculates the risk of neonatal mortality. */
  static void setRiskFromMaternalInfection(int nCounter, int pCounter);
  
  // non-static
  /** Determines whether there is an acute episode or concomitant fever (or
   * neither) and then whether the episode is severe, uncomplicated or there is
   * an indirect death. */
  Pathogenesis::Infection infectionEvent(double ageYears, double totalDensity, double timeStepMaxDensity);
  
  /// Model-specific; for summary.
  virtual double getPyrogenThres();
  
  /// @brief Checkpointing functions
  //@{
  virtual void write(ostream& out) const;
  //@}
  
protected:
  /** Create a PathogenesisModel. */
  PathogenesisModel(double cF);
  /** Create a PathogenesisModel. */
  PathogenesisModel(istream& in);
  
private:	// static
  //comorbidity prevalence at birth as a risk factor for indirect
  static double indirRiskCoFactor_18;
  //sevMal: critical density for severe malaria episode (Y*B1)
  static double sevMal_21;
  //Critical age for co-morbidity (for both severe and indirect)
  static double critAgeComorb_30;
  //comorbidity prevalence at birth as a risk factor for severe
  static double comorbintercept_24;
  
  /** Probability for a newborn to die (indirect death) because the mother is
   * infected. Depends on the prevalence of parasitaemia in mother at some
   * previous t. */
  static double _riskFromMaternalInfection;
  //! array for stored prevalences 20-25 years for 5 months (for neonatal deaths)
  static std::vector<double> _prevalenceByGestationalAge;
  
protected:	// non-static
  virtual double getPEpisode(double timeStepMaxDensity, double totalDensity)=0;
  
  //! comorbidity factor for heterogeneity 
  double _comorbidityFactor; 
};

#endif
