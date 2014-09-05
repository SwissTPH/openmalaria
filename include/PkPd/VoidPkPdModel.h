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

#ifndef Hmod_VoidPkPdModel
#define Hmod_VoidPkPdModel

#include "PkPd/PkPdModel.h"

namespace OM { namespace PkPd {
    /** An implementation of the model which doesn't do anything.
     *
    * Note that there currently needn't be a PKPD model, in which case an instance
    * of this class is created (allowing nicer code). Therefore all methods have an
    * empty implementation in this class.
    */
    class VoidPkPdModel : public PkPdModel {
        virtual void getConcentrations(map<string,double>& concentrations) const{}
        virtual void prescribe(size_t schedule, size_t dosages, double age) {}
        virtual void medicate(double age) {};
        virtual double getDrugFactor (uint32_t proteome_ID) {
            return 1.0;
        }
        virtual void decayDrugs () {}
        virtual uint32_t new_proteome_ID () {
            return 0xFFFFFFFF;
        }
    protected:
        virtual void checkpoint (istream& stream) {}
        virtual void checkpoint (ostream& stream) {}
    };
} }
#endif
