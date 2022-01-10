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

#ifndef Hmod_Infection
#define Hmod_Infection

#include "Global.h"
#include "Parameters.h"
#include "Host/WithinHost/Genotypes.h"

class UnittestUtil;

namespace OM { namespace WithinHost {
    
class Infection {
public:
    inline static void init( SimTime latentP ){
        s_latentP = latentP;
    }
    
    Infection (uint32_t genotype) :
        m_startDate(sim::nowOrTs0()),
        m_density(0.0),
        m_cumulativeExposureJ(0.0),
        m_genotype(genotype)
    {}
    Infection (istream& stream) :
        m_startDate(sim::never())
    {
        m_startDate & stream;
        m_density & stream;
        m_cumulativeExposureJ & stream;
        m_genotype & stream;
    }
    virtual ~Infection () {}
    
    
    /** Return true if infection is blood stage.
     * 
     * Note: infections are considered to be liver stage for 5 days. This is
     * hard-coded since it is convenient in a 5-day time step model (and was one
     * of the reasons a 5-day time step was originally used) (TS).
     * 
     * The remainder of the "latentP" (pre-patent) period is blood-stage, where
     * blood-stage drugs do have an effect but parasites are not detectible.
     * 
     * Note 2: this gets called when deciding which infections to clear. If
     * clearing while updating infections (delayed treatment effect), infections
     * are liver-stage on the time step they start and blood-stage on the next
     * update, thus can be cleared the first time step they are considered
     * blood-stage. If clearing immediately (legacy health system and MDA
     * effect), clearance of blood stage infections can only happen after their
     * first update (though due to the latent period densities will still be
     * low). */
    inline bool bloodStage() const{
        return sim::latestTs0() - m_startDate > sim::fromDays(5);
    }
    
    /// Get the density of the infection as of the last update
    inline double getDensity()const{
        return m_density;
    }
    
    /// Get the cumulative parasite density
    inline double cumulativeExposureJ() {
        return m_cumulativeExposureJ;
    }
    
    /// Get whether the infection is HRP2-deficient
    bool isHrp2Deficient() const {
        return Genotypes::getGenotypes()[m_genotype].hrp2_deficient;
    }
    
    /** Get the infection's genotype. */
    uint32_t genotype()const{ return m_genotype; }
    
    
    /// Resets immunity properties specific to the infection (should only be
    /// called along with clearImmunity() on within-host model).
    inline void clearImmunity(){
        m_cumulativeExposureJ = 0.0;
    }
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        checkpoint (stream);
    }
    
protected:
    inline virtual void checkpoint (ostream& stream) {
        m_startDate & stream;
        m_density & stream;
        m_cumulativeExposureJ & stream;
        m_genotype & stream;
    }
    
    /// Date of inoculation of infection (start of liver stage)
    /// This is the step of inoculation (ts0()).
    SimTime m_startDate = sim::never();
    
    /// Current density of the infection
    double m_density;
    
    /// Cumulative parasite density, since start of this infection
    double m_cumulativeExposureJ;
    
private:
    /// Genotype of infection (a code; see Genotypes class).
    uint32_t m_genotype;
    
protected:
    /// Pre-erythrocytic latent period (instantiated in WHFalciparum.cpp)
    static SimTime s_latentP;
    
    friend class ::UnittestUtil;
};

} }
#endif
