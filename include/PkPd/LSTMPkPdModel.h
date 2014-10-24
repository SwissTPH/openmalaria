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

namespace scnXml{
    class PKPDMedication;
}
namespace OM { namespace PkPd {

struct MedicateData {
    MedicateData () :
        drug(0),
        qty(numeric_limits< double >::signaling_NaN()),
        time(numeric_limits< double >::signaling_NaN()),
        duration(numeric_limits< double >::quiet_NaN())
    {}
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        drug & stream;
        qty & stream;
        time & stream;
        duration & stream;
    }
    
private:
    void load( const scnXml::PKPDMedication& med );
    
    inline MedicateData multiplied( double doseMult ){
        MedicateData r( *this );
        r.qty *= doseMult;
        return r;
    }
    
    size_t drug;      /// Drug type index
    double qty;         /// Quantity of drug prescribed (mg when oral, mg/kg when IV)
    double time;        /// Time to medicate at, in days (0 means start of time step, may be >= 1 (thus not today))
    double duration;    /// Duration for IV purposes, in days (use 0 or NaN to indicate oral dose)
    
    friend class Schedule;
    friend class LSTMPkPdModel;
    friend class ::UnittestUtil;
};

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
    
    virtual void getConcentrations(map<string,double>& concentrations) const;
    virtual void prescribe(size_t schedule, size_t dosages, double age, double body_mass);
    virtual void medicate(double body_mass);
    virtual double getDrugFactor (uint32_t proteome_ID);
    virtual void decayDrugs ();
    
    virtual uint32_t new_proteome_ID ();
    
private:
  /** Medicate drugs to an individual, which act on infections the following
   * time steps, until rendered ineffective by decayDrugs().
   *
   * \param typeIndex The index of drug type data (what LSTMDrugType::findDrug() returns).
   * \param qty The quantity in either mg (if oral dose) or mg/kg (if IV).
   * \param time Time in days since start of this time step to medicate at
   * \param duration  Duration in days. 0 or an NaN indicates no duration.
   * \param bodyMass Weight of human in kg
   * 
   * Due to the fact we're using a discrete time step model, the case-management
   * update (calling medicate) and within-host model update (calling
   * getDrugFactor) cannot [easily] have immediate effects on each other. The
   * implementation we use is that the within-host model update (calculating
   * new infection densities) happens first; hence medicate() will always be
   * called after getDrugFactor in a time step, and a time of zero means the
   * dose has effect from the start of the following time step. */
  void medicateDrug(size_t typeIndex, double qty, double time, double duration, double bodyMass);
    
    /// Drugs with non-zero blood concentrations:
    list<LSTMDrug> _drugs;
    
    /// All pending medications
    list<MedicateData> medicateQueue;
    
    friend class ::UnittestUtil;
};

} }
#endif