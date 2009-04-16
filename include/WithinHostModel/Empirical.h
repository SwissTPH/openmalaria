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

class EmpiricalWithinHostModel : public Infection {
public:
  EmpiricalWithinHostModel ();
  double getNewDensity(double * transformedLaggedDensities, int ageOfInfection);
  void initialiseInfection(double * transfomedLaggedDensities);
  void setInflationFactors(double inflationMean, double inflationVariance);
private:
  double getInflatedDensity(double nonInflatedDensity);
  double sigma_noise(int ageOfInfection);
  double samplePatentValue(double mu, double sigma, double lowerBound);
  double sampleSubPatentValue(double mu, double sigma, double upperBound);
  double inverseBoxCoxTransform(double transformedValue);
  double boxCoxTransform(double untransformedValue);
  static const int _maximumDurationInDays=418; 
  double _maximumPermittedAmplificationPerCycle;
  double _subPatentLimit;
  double _inflationVariance;
  double _inflationMean;
  double _lambda;
  double _alpha1;
	double _alpha2;	
  double _alpha3;
  double _sigma_alpha1;	
	double _sigma_alpha2;	
	double _sigma_alpha3;
  double _sigma0_res;	
  double _sigmat_res;
  double _mu_beta1[_maximumDurationInDays];
  double _sigma_beta1[_maximumDurationInDays];
  double _mu_beta2[_maximumDurationInDays];
  double _sigma_beta2[_maximumDurationInDays];
  double _mu_beta3[_maximumDurationInDays];
  double _sigma_beta3[_maximumDurationInDays];  
};
