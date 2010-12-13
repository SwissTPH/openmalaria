/*
 This file is part of OpenMalaria.

 Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine

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

#ifndef Hmod_AgeDistribution
#define Hmod_AgeDistribution

#include "Global.h"
#include "inputData.h"

#include <limits>
#include <stdexcept>

namespace OM { namespace util {
    
/** A class representing deterministic interpolation of data collected
 * according to age groups. Derived classes implement the actual interpolation.
 * 
 * Current version does not store an age index, thus an order log(n) lookup
 * must occur each time a value is looked up.
 ********************************************/
class AgeGroupInterpolation
{
public:
    /** Return a dummy object (avoids dangling pointer). */
    static AgeGroupInterpolation* dummyObject();
    /** Return a new age-group-data interpolator. XML fragment specifies
     * which interpolation and values to use.
     * 
     * @param ageGroups XML element of per-age-group values
     * @param eltName Name of XML element (to all reasonable error reporting)
     */
    static AgeGroupInterpolation* makeObject(
        const scnXml::AgeGroupValues& ageGroups, const char* eltName
    );
    /** Free pointed object.
     *
     * (Note: dummy object should not be freed; this checks for that.) */
    static void freeObject( AgeGroupInterpolation* obj );
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        checkpoint( stream );
    }
    
    /** Return a value interpolated for age ageYears. */
    virtual double operator() (double ageYears) const =0;
    
    /** Scale held fatality rate by factor. */
    virtual void scale( double factor ) =0;
    
protected:
    virtual void checkpoint (ostream& stream) =0;
    virtual void checkpoint (istream& stream) =0;
};


/** This class gives direct access to input age-group
 * data (discontinuous).
 ********************************************/
class AgeGroupPiecewiseConstant : public AgeGroupInterpolation
{
public:
    AgeGroupPiecewiseConstant (
        const scnXml::AgeGroupValues& ageGroups, const char* eltName
    );
    
    virtual double operator() (double ageYears) const;
    
    virtual void scale( double factor );
    
protected:
    virtual void checkpoint (ostream& stream);
    virtual void checkpoint (istream& stream);
    
    // Points to interpolate between in the middle of input age groups. Extra
    // points at zero and infinity are added with equal value to first and last
    // points respectively.
    map<double,double> dataPoints;
};


/** This class gives piecewise linear inpolation on top of input age-group
 * data (continuous but with discontinuous derivative).
 ********************************************/
class AgeGroupPiecewiseLinear : public AgeGroupInterpolation
{
public:
    AgeGroupPiecewiseLinear (
        const scnXml::AgeGroupValues& ageGroups, const char* eltName
    );
    
    virtual double operator() (double ageYears) const;
    
    virtual void scale( double factor );
    
protected:
    virtual void checkpoint (ostream& stream);
    virtual void checkpoint (istream& stream);
    
    // Points to interpolate between in the middle of input age groups. Extra
    // points at zero and infinity are added with equal value to first and last
    // points respectively.
    map<double,double> dataPoints;
};


} }
#endif