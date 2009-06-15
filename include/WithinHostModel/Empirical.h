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

#include "WithinHostModel/Infection.h"

class EmpiricalInfection : public Infection {
public:
  EmpiricalInfection(double growthRateMultiplier);
  static void initParameters();
  double getNewDensity(int ageOfInfection, double growthRateMultiplier);
  static void overrideInflationFactors(double inflationMean, double inflationVariance, double extinctionLevel, double overallMultiplier);


private:
  double getInflatedDensity(double nonInflatedDensity);
  double sigma_noise(int ageOfInfection);
  double samplePatentValue(double mu, double sigma, double lowerBound);
  double sampleSubPatentValue(double mu, double sigma, double upperBound);
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
  double _laggedLogDensities[3];
  int _startTime;
};
