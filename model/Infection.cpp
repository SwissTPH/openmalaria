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

#include "Infection.h"

float Infection::cumulativeYstar;
float Infection::cumulativeHstar;

double sampleFromLogNormal(double normp, double meanlog, double stdlog){
    /*
    
  Used for performance reasons. Calling GLS LOG_NORMAL 5 times is 50% slower.
    
    */

  double zval;
  double valsampleFromLogNormal;
  zval=W_UGAUSS_PINV(normp);
    /*
  Why not zval=W_UGAUSS?
  where normp is distributed uniformly over [0,1],
  zval is distributed like a standard normal distribution
  where normp has been transformed by raising to the power of 1/(T-1) 
  zval is distributed like a uniform gauss	times 4* F(x,0,1)^3, where F(x,0,1) ist the cummulative
  distr. function of a uniform gauss
    */
  valsampleFromLogNormal=exp(meanlog+stdlog*((float)zval));
  return valsampleFromLogNormal;
}
