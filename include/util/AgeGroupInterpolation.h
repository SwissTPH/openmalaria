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

namespace scnXml {
    class AgeGroupValues;
}
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
    virtual ~AgeGroupInterpolation() {}
    
    /** Return a dummy object (avoids dangling pointer). */
    static AgeGroupInterpolation* dummyObject();
    /** Return a new age-group-data interpolator. XML fragment specifies
     * which interpolation and values to use.
     * 
     * @param ageGroups XML element of per-age-group values
     * @param eltName Name of XML element (for reasonable error reporting)
     */
    static AgeGroupInterpolation* makeObject(
        const scnXml::AgeGroupValues& ageGroups, const char* eltName
    );
    /** Free pointed object.
     *
     * (Note: dummy object should not be freed; this checks for that.) */
    static void freeObject( AgeGroupInterpolation* obj );
    
    /// Return true if instance represents something other than the dummy object.
    inline bool isSet() {
        return this != dummyObject();
    }
    
    /** Return a value interpolated for age ageYears. */
    virtual double eval (double ageYears) const =0;
    
    /** Scale function by factor. */
    virtual void scale( double factor ) =0;

    /** Find the youngest age which is the global maximum (i.e. the age at
     * which individuals are considered adults where all adults have the same
     * same value). */
    virtual double firstGlobalMaximum() =0;
    
protected:
    /** Sample interpolator between 0 and max age, outputting to a csv file
     * called name.csv. */
    void outputSamples( const string name );
};

} }
#endif
