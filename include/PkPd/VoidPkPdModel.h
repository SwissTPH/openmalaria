/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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
        virtual void medicate(string drugAbbrev, double qty, double time, double ageYears) {}
        virtual void medicateIV(string drugAbbrev, double qty, double duration, double endTime) {}
        virtual void decayDrugs () {}
        virtual double getDrugFactor (uint32_t proteome_ID) {
            return 1.0;
        }
        virtual uint32_t new_proteome_ID () {
            return 0xFFFFFFFF;
        }
    protected:
        virtual void checkpoint (istream& stream) {}
        virtual void checkpoint (ostream& stream) {}
    };
} }
#endif
