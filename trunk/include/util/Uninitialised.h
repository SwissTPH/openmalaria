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

#ifndef UTIL_UNINITIALISED
#define UTIL_UNINITIALISED

#include "util/errors.h"

namespace OM { namespace util {

/** Wrapper which throws if variables are used before initialisation.
 * 
 * Usage: replace "T var" with "Uninitialised<T> var". If "var"'s value is used
 * before it is set, the wrapper will throw. Note: variables containing
 * Not-A-Number values are considered uninitialised.
 * 
 * This is intended for simple types only (only assignment and cast supported).
 * Only intended for debugging.
 */
template<typename T>
class Uninitialised {
public:
    Uninitialised () : initialised(false) {}
    Uninitialised (T var) : variable(var) {
        // don't consider setting to NaN initialisation
        initialised = (var == var);
    }
    
    inline void operator= (T var) {
        variable = var;
        initialised = (var == var);
    }
    
    // Implicit cast-to-T (should be invoked automatically on use): throw if uninitialised
    inline operator T () const{
        if( !initialised )
            throw traced_exception("uninitialised variable used!");
        return variable;
    }
    
private:
    T variable;
    bool initialised;
};

} }

#endif

