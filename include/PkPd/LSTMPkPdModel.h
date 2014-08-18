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

#ifndef Hmod_LSTMPkPdModel
#define Hmod_LSTMPkPdModel

#include "PkPd/PkPdModel.h"
#include "PkPd/Drug/LSTMDrug.h"
#include "PkPd/LSTMTreatments.h"

namespace OM { namespace PkPd {
    
/** Pharmacokinetic and pharmacodynamics interface, used by each human's
 * within-host model.
 *
 * Some of the implementation is contained in the drug.h/drug.cpp files.
 * 
 * This class holds per-human data: prescribed medications and drugs in the
 * body. */
class LSTMPkPdModel : public PkPdModel {
public:
    LSTMPkPdModel ();
    virtual ~LSTMPkPdModel ();
    
    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);
    
    virtual void prescribe(size_t schedule, size_t dosages, double age);
    
    //TODO: call
    /// Call if there are any pending medications.
    void doUpdate( double bodyMass );
    
    virtual void decayDrugs ();
    virtual double getDrugFactor (uint32_t proteome_ID);
    
    virtual uint32_t new_proteome_ID ();
    
private:
  /** Medicate drugs to an individual, which act on infections the following
   * timesteps, until rendered ineffective by decayDrugs().
   *
   * \param typeIndex The index of drug type data (what LSTMDrugType::findDrug() returns).
   * \param qty The quantity in mg.
   * \param time Time in days since start of this time step to medicate at
   * \param duration  Duration in days. 0 or an NaN indicates no duration.
   * \param bodyMass Weight of human in kg
   * 
   * Due to the fact we're using a discrete timestep model, the case-management
   * update (calling medicate) and within-host model update (calling
   * getDrugFactor) cannot [easily] have immediate effects on each other. The
   * implementation we use is that the within-host model update (calculating
   * new infection densities) happens first; hence medicate() will always be
   * called after getDrugFactor in a timestep, and a time of zero means the
   * dose has effect from the start of the following timestep. */
  void medicate(size_t typeIndex, double qty, double time, double duration, double bodyMass);
    
    /// Drugs with non-zero blood concentrations:
    list<LSTMDrug> _drugs;
    
    /// All pending medications
    list<MedicateData> medicateQueue;
    
    friend class ::UnittestUtil;
};

} }
#endif