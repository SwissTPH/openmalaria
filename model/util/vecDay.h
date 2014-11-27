/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2014 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2014 Liverpool School Of Tropical Medicine
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

#ifndef Hmod_OM_util_vecDay
#define Hmod_OM_util_vecDay

#include "Global.h"
#include <vector>

#if __cplusplus >= 201103L
#error "vecDay should be updated if switching to a later C++ version"
#endif

namespace OM {
namespace util {

/** A std::vector<T,Alloc>, but where the index is in days (via SimTime type).
 * 
 * Rather hacky, with only the needed bits implemented. */
template<typename T, typename Alloc = std::allocator<T> >
struct vecDay {
    typedef std::vector<T, Alloc> vec_t;
    
    vecDay() : v() {}
    explicit vecDay(const typename vec_t::allocator_type& a) : v(a) {}
    explicit vecDay(SimTime n, const typename vec_t::value_type& value =
        typename vec_t::value_type(), const typename vec_t::allocator_type& a =
        typename vec_t::allocator_type() ) :
        v(static_cast<size_t>(n.inDays()), value, a) {}
    vecDay(const vecDay& x) : v(x.v) {}
    
    inline void assign(SimTime n, const typename vec_t::value_type& val){
        v.assign(n.inDays(), val); }

    inline void resize(SimTime new_size, typename vec_t::value_type x = typename vec_t::value_type()){
        v.resize( new_size.inDays(), x ); }
    
    inline SimTime size() const {
        return sim::fromDays(v.size()); }
    
    inline typename vec_t::reference
    operator[](SimTime n){ return v[n.inDays()]; }
    
    inline typename vec_t::const_reference
    operator[](SimTime n) const{ return v[n.inDays()]; }
    
    /// Access
    const vec_t& internal()const{ return v; }
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        v & stream;
    }
    
private:
    vec_t v;
};

/** A two-dimensional array whose first index has SimTime type (as vecDay) and
 * second index has size_t type. */
template<typename T, typename Alloc = std::allocator<T> >
struct vecDay2D {
    typedef std::vector<T, Alloc> vec_t;
    
    vecDay2D() : v() {}
    explicit vecDay2D(const typename vec_t::allocator_type& a) : stride(0), v(a) {}
    explicit vecDay2D(SimTime n1, size_t n2,
        const typename vec_t::value_type& value = typename vec_t::value_type(),
        const typename vec_t::allocator_type& a = typename vec_t::allocator_type() )
            : stride(n2), v(static_cast<size_t>(n1.inDays() * n2), value, a) {}
    vecDay2D(const vecDay2D& x) : stride(x.stride), v(x.v) {}
    
    inline void assign(SimTime dim1, size_t dim2,
        const typename vec_t::value_type& val)
    {
        v.assign(dim1.inDays() * dim2, val);
        stride = dim2;
    }

    inline void resize(SimTime dim1, size_t dim2,
        typename vec_t::value_type x = typename vec_t::value_type())
    {
        v.resize( dim1.inDays() * dim2, x );
        stride = dim2;
    }
    
    inline typename vec_t::reference
    at(SimTime n1, size_t n2){
        return v[n1.inDays() * stride + n2];
    }
    
    inline typename vec_t::const_reference
    at(SimTime n1, size_t n2) const{
        return v[n1.inDays() * stride + n2];
    }
    
    inline vec_t& internal_vec(){ return v; }
    
    inline void set_all( typename vec_t::value_type x ){
        v.assign( v.size(), x );
    }
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        stride & stream;
        v & stream;
    }
    
private:
    size_t stride;
    vec_t v;
};

namespace vectors{
  /// Scale all elements of a vector by a in-situ
  void scale (vecDay<double>& vec, double a);
  
  /// Scale all elements of a vector by a in-situ
  void scale (vecDay2D<double>& vec, double a);
  
  /// Return sum of all elements
  double sum (const vecDay<double>& vec);
  
  /** The inverse of logDFT (or an approximation, when N&lt;T or
   * tArray.size() ≠ T). Result may also be rotated.
   * 
   * (This was called calcExpFourierSeries, and does essentially the same
   * thing.)
   *
   * @param tArray Array to fill with exponated values from Fourier series.
   * Length should already be set. Need not have the same length as the
   * array used to calculate FC.
   * @param FC Fourier coefficients (a0, a1,b1, a2,b2, ...); can be any length
   * so long as it is odd.
   * @param rAngle Angle to rotate generated series by in radians: [0,2π] */
  void expIDFT (vecDay<double>& tArray, const vector<double>& FC, double rAngle);
}

/// Utility to print a vector (operator must be in namespace)
template<class T>
ostream& operator<< (ostream& out, vecDay<T> vec) {
    out << '[';
    if (vec.size() > sim::zero())
        out << vec[sim::zero()];
    for( SimTime i = sim::oneDay(), end = vec.size(); i < end; i += sim::oneDay() ){
        out << ", " << vec[i];
    }
    out << ']';
    return out;
}

}
}
#endif
