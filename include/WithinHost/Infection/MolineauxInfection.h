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

#ifndef Hmod_MOLINEAUXINFECTION_H
#define Hmod_MOLINEAUXINFECTION_H

#include "WithinHost/Infection/CommonInfection.h"

class MolineauxInfectionSuite;

namespace OM { namespace WithinHost {

/**
 * Implementation of a P. falciparum model by Molineaux et. al., published in:
 * 
 * L. MOLINEAUX, H. H. DIEBNER, M. EICHNER, W. E. COLLINS, G. M. JEFFERY and
 * K. DIETZ, 2001: Plasmodium falciparum parasitaemia described by a new
 * mathematical model. Parasitology, 122, pp 379-391
 * doi:10.1017/S0031182001007533
 */
class MolineauxInfection : public CommonInfection {
public:
    ///@brief Static class members
    //@{
    static void init(const OM::Parameters& parameters);
    
    // Constants:
    // v: number of variants per clone (one infection = one new clone); 50 in paper
    //TODO: v should have a significant effect on performance, but reducing by,
    // say, half, may not have a big effect on the model. Evaluate.
    static const size_t v = 50;
    // taus: used for the variantTranscending and variantSpecific array, 4 Molineaux time steps = 8 days
    static const size_t taus = 4;
    //@}
    
    // Initialise. Samples several parameters.
    MolineauxInfection(uint32_t protID);
    // Load from a checkpoint:
    MolineauxInfection (istream& stream);
    virtual ~MolineauxInfection () {};
    
    virtual bool updateDensity( double survivalFactor, SimTime bsAge, double body_mass );
    
protected:
    virtual void checkpoint (ostream& stream);
    
private:
    double getVariantSpecificSummation(int i, double P_current);
    
    // NOTE: we also have inherited parameters
    // m_startDate is used to give the age here
    // m_density is equivalent to Pc in paper
    // m_cumulativeExposureJ is cumulative parasite density (used by external immunity function)
    
    // m[i]: Multiplication factor, per two-day cycle of variant i
    float m[v];
    // variantTranscendingSummation: See Molineaux paper, equation 7
    float variantTranscendingSummation;
    // index: we use (ageDays mod 8) / 2 (but ageDays could be replaced by e.g. simTimeDays if available)
    float laggedPc[taus];
    /* Pc_star, Pm_star: two host-specific critical densities.
     * These two values depend on the first local maximum or the difference
     * between the last positive day and the first positive day. */
    float Pc_star, Pm_star;
    
    struct Variant {
        // Initialise all variables to 0
        Variant ();
        
        /// Checkpointing
        void operator& (ostream& stream);
        void operator& (istream& stream);
        
        //TODO: we probably don't need P now
        // P_i(t+1), P_i(t+2): variant's i density (PRBC/μl blood)
        float P1, P2;
        // variantSpecificSummation: See Molineaux paper, equation 6
        float variantSpecificSummation;
        // index: we use (ageDays mod 8) / 2 (but ageDays could be replaced by e.g. simTimeDays if available)
        float laggedP[taus];
    };
    // variant-specific data; variants[i-1] corresponds to variant i in the paper
    //TODO: more optimal storage of variants? Fixed-size list? Like this, or a map?
    vector<Variant> variants;
    
    // allow unittest to access private vars
    friend class ::MolineauxInfectionSuite;
};

}}
#endif
