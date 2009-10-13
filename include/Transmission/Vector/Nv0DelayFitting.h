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

#ifndef Hmod_VectorNv0DelayFitting
#define Hmod_VectorNv0DelayFitting

#include <boost/foreach.hpp>
#include <boost/math/tools/roots.hpp>

namespace Nv0DelayFitting {
template <class T>
struct eDFunctor
{
  eDFunctor(double r, const vector<T>& fc_, const vector<T>& samples_) : p(samples_.size()), existingR(r), fc(fc_), logSamples(samples_) {
    if (fc.size() % 2 == 0)
      throw runtime_error("The number of Fourier coefficents should be odd.");
    w = 2*M_PI / T(p);
    fn = (fc.size()-1)/2;
    BOOST_FOREACH (T& sample, logSamples) {
      sample = log (sample);	// compare logarithms of EIR to make differentiation easier
    }
  }
  
  std::tr1::tuple<T, T, T> operator()(T const& d)
  {
    T f = 0.0, df = 0.0, ddf = 0.0;
    
    // Calculate inverse discrete Fourier transform
    for (size_t t=0; t<p; t++) {
      T wt = w*t+d;
      T val = fc[0], dval = 0.0, ddval = 0.0;
      for(size_t n=1;n<=fn;n++){
	T temp = fc[2*n-1]*cos(n*wt) + fc[2*n]*sin(n*wt);	// value
	val  += temp;
	dval += n*fc[2*n]*cos(n*wt) - n*fc[2*n-1]*sin(n*wt);	// first derivative wrt d
	ddval = -n*n*temp;					// 2nd derivative wrt d
      }
      
      // The difference of logarithms of sample and fourier value
      T diff = val - logSamples[t];		// deriv. wrt. d is just dval
      f += diff*diff;				// add diff²
      df += 2*diff * (dval - 0.0);		// add 1st deriv. diff²
      ddf += 2*dval*dval + 2*diff*ddval;	// add 2nd deriv. diff²
      cout << t << "\t"<< val << "\t" << logSamples[t] << "\t" << diff*diff << "\t" <<2*diff * (dval - 0.0) << "\t" << 2*dval*dval + 2*diff*ddval << endl;
    }
    
    return std::tr1::make_tuple(f, df, ddf);
  }
  private:
    size_t p, fn;
    T w,existingR;
    const vector<T>& fc;
    vector<T> logSamples;
};


/** Calculate rotation angle needed to match up fourier series defined by fc with samples.
 *
 * @param existingR Existing angle (in radians) to rotate by.
 * @param fc Fourier coefficients for S_v series (a0, a1, b1, ...)
 * @param samples The calculated S_v values we want to match
 * @returns Angle to rotate by (including existingR).
 */
template <class T>
T fit (double existingR, const vector<T>& fc, const vector<T>& samples)
{
  T min = 0.0;
  T max = 2*M_PI;	// can we provide better bounds?
  T guess = 0.0;
  int digits = std::numeric_limits<T>::digits / 2;	//TODO: check how many digits we should provide
  //TODO: see whether halley_iterate or schroeder_iterate is faster
  return -boost::math::tools::halley_iterate(eDFunctor<T>(existingR, fc, samples), guess, min, max, digits) + existingR;
}
}
#endif
