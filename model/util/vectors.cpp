/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#include "util/vectors.h"
#include "util/errors.h"
#include <cstring>
#include <cmath>
#include <cassert>

namespace OM { namespace util {
    
void vectors::scale (vector<double>& vec, double a) {
  for (size_t i = 0; i < vec.size(); ++i)
    vec[i] *= a;
}

double vectors::sum (const vector<double>& vec) {
  double r = 0.0;
  for (size_t i = 0; i < vec.size(); ++i)
    r += vec[i];
  return r;
}

double vectors::sum (const gsl_vector *vec) {
  double r = 0.0;
  for (size_t i = 0; i < vec->size; ++i)
    r += gsl_vector_get( vec, i );
  return r;
}

void vectors::addTo (vector<double>& x, vector<double>& y){
    assert( x.size() == y.size() );
    for( size_t i=0; i<x.size(); ++i ){
	x[i] += y[i];
    }
}


bool vectors::approxEqual (const double a, const double b, const double lim_fact) {
    bool aE = (fabs(a-b) <= max(fabs(a),fabs(b)) * lim_fact);
#ifndef NDEBUG
    if (!aE){
        cerr<<"not approx equal: "<<a<<", "<<b<<endl;
    }
#endif
    return aE;
}

bool vectors::approxEqual (const vector<double>& vec1, const vector<double>& vec2, const double lim_fact) {
  if (vec1.size() != vec2.size())
    return false;
  for (size_t i = 0; i < vec1.size(); ++i) {
    if (!approxEqual (vec1[i], vec2[i], lim_fact))
      return false;
  }
  return true;
}


vector<double> vectors::gsl2std (const gsl_vector* vec) {
  vector<double> ret (vec->size);
  memcpy (&ret[0], vec->data, ret.size()*sizeof(double));
  return ret;
}
void vectors::gsl2std( const gsl_vector *vec, vector<double>& target ){
  target.resize (vec->size);
  memcpy (&target[0], vec->data, target.size()*sizeof(double));
}


gsl_vector* vectors::std2gsl (const vector<double>& vec, size_t length) {
  if (vec.size() != length)
    throw TRACED_EXCEPTION_DEFAULT ("vectorStd2Gsl: vec has incorrect length");
  gsl_vector* ret = gsl_vector_alloc (length);
  memcpy (ret->data, &vec[0], length * sizeof(double));
  return ret;
}

gsl_vector* vectors::std2gsl (const double* vec, size_t length) {
  gsl_vector* ret = gsl_vector_alloc (length);
  memcpy (ret->data, vec, length * sizeof(double));
  return ret;
}

void vectors::logFourierCoefficients(const vector<double>& iArray, vector<double>& FC) {
    logDFT(iArray, FC);
    for (size_t i=1; i<FC.size(); ++i)
        FC[i] *= 2;
}

// NOTE: could replace these algorithms with FFTs, but since these aren't
// called in performance-critical code, the difference would be insignificant.
void vectors::logDFT(const vector<double>& iArray, vector<double>& FC) {
    if (FC.size() > 2*iArray.size() || mod_nn(FC.size(), 2) == 0)
        throw TRACED_EXCEPTION_DEFAULT("Require DFT FC.size() to be 2*n-1 for n<=iArray.size()");
    
    size_t T = iArray.size();
    size_t N = (FC.size() + 1) / 2;
    
    double w = 2.0 * M_PI / T;
    FC.assign(FC.size(), 0.0);
    
    for (size_t t = 0; t < T; ++t) {
        double val = log( iArray[t] );
        FC[0] += val;
        for (size_t n = 1; n < N; ++n ){
            FC[2*n-1] += val * cos(w*n*t);
            FC[2*n  ] += val * sin(w*n*t);
        }
    }
    scale(FC, 1.0 / T);
}

void vectors::expIDFT (std::vector< double >& tArray, const std::vector< double >& FC, double rAngle) {
    if (mod_nn(FC.size(), 2) == 0)
        throw TRACED_EXCEPTION_DEFAULT("The number of Fourier coefficents should be odd.");
    
    size_t T2 = tArray.size();
    size_t N = (FC.size() + 1) / 2;
    
    double w = 2.0 * M_PI / T2;
    
    // Calculate inverse discrete Fourier transform
    // TODO: This may not interpolate sensibly. See for example
    // https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Trigonometric_interpolation_polynomial
    for (size_t t = 0; t < T2; ++t) {
        double temp = FC[0];
        double wt = w*t - rAngle;
        for (size_t n = 1; n < N; ++n) {
            temp += FC[2*n-1]*cos(n*wt) + FC[2*n]*sin(n*wt);
        }
        tArray[t] = exp(temp);
    }
}


} }