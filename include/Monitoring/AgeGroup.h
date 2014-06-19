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

#ifndef Hmod_Monitoring_AgeGroup
#define Hmod_Monitoring_AgeGroup

#include "Global.h"
#include "util/errors.h"

namespace scnXml{ class Monitoring; }
namespace OM {
namespace Monitoring {
/**
 * Included for type-saftey: don't allow implicit double->int conversions.
 *
 * Incindentally, the constructor can be used implicitly for implicit
 * conversion doing the right thing.
 * 
 * Don't use _this_ class for other index/age-group types. */
class AgeGroup {
  public:
    AgeGroup () : index(0) {}
    
    /** Update age-group. Assumes age only increases (per instance).
     *
     * If called regularly, should be O(1); worst case is O(_upperbound.size()). */
    void update (double ageYears);
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        index & stream;
    }
    
    /** Get the represented index. */
    inline size_t i () const{
        return index;
    }
    
    /// Get the total number of age categories (inc. one for indivs. not in any
    /// category given in XML).
    static inline size_t getNumGroups () {
        if( _upperbound.size() == 0 ) throw TRACED_EXCEPTION_DEFAULT( "not yet initialised" );
        return _upperbound.size();
    }
    
private:
    size_t index;
    
    /// Read age group bounds from XML data
    static void init (const scnXml::Monitoring& monitoring);
    
    //BEGIN Static parameters only set by init()
    /** Upper boundary of agegroups, in years.
     *
     * These are age-groups given in XML plus one with no upper limit for
     * individuals outside other bounds. */
    static vector<double> _upperbound;
    //END
    
    friend class Survey;
};

} }
#endif
