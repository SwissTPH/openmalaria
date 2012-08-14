/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2011 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#ifndef HS_DIAGNOSTIC
#define HS_DIAGNOSTIC

#include "Global.h"
#include "schema/interventions.h"
#include <limits>

namespace OM { namespace Clinical {
    
    class Diagnostic {
    public:
        /** Construct. Set params to NaNs. */
        Diagnostic() :
            specificity( numeric_limits<double>::signaling_NaN() ),
            density( numeric_limits<double>::signaling_NaN() )
        {}
        
        /** Set parameters from an XML element. */
        void init( const scnXml::HSDiagnostic& elt );
        
        /** Use the test.
         * 
         * @param x Current parasite density in parasites per µL
         * @returns True if outcome is positive. */
        bool isPositive( double x ) const;
        
    private:
        // either not-a-number indicating a deterministic test, or specificity
        double specificity;
        // depending on model, this is either the minimum detectible density
        // or the density at which test has half a chance of a positive outcome
        double density;
    };

} }
#endif
