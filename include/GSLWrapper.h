/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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


#ifndef GSL_WRAPPER_H
#define GSL_WRAPPER_H

double W_BETA(double a, double b);
double W_GAUSS(double mean, double std);
double w_minimize_calc_rss(double* par1, double* par2);
double W_UGAUSS_P(double x);
double W_UGAUSS_PINV(double p);
double W_LOGNORMAL(double mean, double std);
int W_POISSON(double lambda);
double W_GAMMA(double a, double b);
void save_rng_state(int seedFileNumber);
void load_rng_state(int seedFileNumber);
double W_UNIFORM();
double SETDEMOPARAMETERS(double *param1,double*param2);

void  GSL_SETUP();
void  GSL_TEARDOWN();

#endif
