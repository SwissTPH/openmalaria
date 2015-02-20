/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2015 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
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

#ifndef Hmod_LSTMDrugThreeComp
#define Hmod_LSTMDrugThreeComp

#include "Global.h"
#include "PkPd/Drug/LSTMDrug.h"
#include "PkPd/Drug/LSTMDrugType.h"

using namespace std;

namespace OM {

namespace PkPd {
struct DoseParams;
}
namespace PkPd {

/** A class holding PK and PD drug info, per human, per drug type.
 * 
 * Three compartment model:
 * 
 * "Mathematical Expressions of the Pharmacokinetic and Pharmacodynamic Models
 * implemented in the Monolix software", Julie Bertrand and France Mentré, Sept
 * 2008. From page 37, in particular equation (1.72) p44. */
class LSTMDrugThreeComp : public LSTMDrug {
public:
    /** Create a new instance. */
    LSTMDrugThreeComp (const LSTMDrugType&);
    
    virtual size_t getIndex() const;
    virtual double getConcentration() const;
    
    virtual void medicate (double time, double qty, double bodyMass);
    
    virtual double calculateDrugFactor(uint32_t genotype) const;
    virtual void updateConcentration ();
    
protected:
    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);
    
    /// Always links a drug instance to its drug-type data
    //TODO: does it still make sense to link this?
    const LSTMDrugType& typeData;
    
    // Concentrations in the blood and other compartments
    // See "Permutation of the three-compartment equation",
    // Diggory Hardy, Swiss TPH, February 3, 2015
    //NOTE: We assume instantaneous absorbtion, thus do not need the ABC term
    // (rationale: lets not complicate code without evidence of need)
    double concA, concB, concC, concABC;
    
    // Computed sampled constants
    double na, nb, ng, nka;      // -α, -β, -γ, -k_a
    double A, B, C;
};

}
}
#endif
