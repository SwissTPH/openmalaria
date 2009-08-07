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

#include "WithinHost/Infection.h"

class EmpiricalInfection : public Infection {
public:
  /// Per instance initialisation; create new inf.
  EmpiricalInfection(int startTime, double growthRateMultiplier);

  /** Checkpoint-reading constructor */
  EmpiricalInfection (istream& in);

  /** Destructor
   * 
   * Note: this destructor does nothing in order to allow shallow copying to
   * the population list. destroy() does the real freeing and must be
   * called explicitly. */
  ~EmpiricalInfection();
  void destroy();


  void write (ostream& out) const;

  static void init ();

  //! Get the last timestep before the infection is cleared.
  /*!
    \return The interval before clearance.
  */
  int getEndDate();

  /// Multiplies the density by x.
  void multiplyDensity(double x) { _density *= x; };
  //! Get the density of the infection
  double getDensity() { return _density; };

  //! Start date of the infection
  int getStartDate() { return _startdate; };


  //! Set patent growth rate multiplier
  void setPatentGrowthRateMultiplier(double multiplier);

  int getDuration() { return _duration; };



  /// Static (shared) data initialisation
  static void initParameters();
  
  /// update
  double determineWithinHostDensity();
  
  /// only for parameterisation?
  static void overrideInflationFactors(double inflationMean, double inflationVariance, double extinctionLevel, double overallMultiplier);


private:
  double getInflatedDensity(double nonInflatedDensity);
  double sigma_noise(int ageOfInfection);
  double samplePatentValue(double mu, double sigma, double lowerBound);
  double sampleSubPatentValue(double mu, double sigma, double upperBound);
  //the following were meant to be static but this lead to a linker error
  static const int _maximumDurationInDays=418; 
  static double _maximumPermittedAmplificationPerCycle;
  static double _subPatentLimit;
  static double _lambda;
  static double _alpha1;
  static double _alpha2;	
  static double _alpha3;
  static double _mu1;	
  static double _mu2;	
  static double _mu3;
  static double _sigma0_res;	
  static double _sigmat_res;
  static double _mu_beta1[_maximumDurationInDays];
  static double _sigma_beta1[_maximumDurationInDays];
  static double _mu_beta2[_maximumDurationInDays];
  static double _sigma_beta2[_maximumDurationInDays];
  static double _mu_beta3[_maximumDurationInDays];
  static double _sigma_beta3[_maximumDurationInDays];
  static double _inflationMean;
  static double _inflationVariance;
  static double _extinctionLevel;
  static double _overallMultiplier;
  // the following are not static
  double _laggedLogDensities[3];
  int _startTime;
  double _patentGrowthRateMultiplier;
};
