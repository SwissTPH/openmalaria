/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2013 Swiss Tropical and Public Health Institute 
 * Copyright (C) 2005-2013 Liverpool School Of Tropical Medicine
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

#ifndef OM_INTERVENTIONS_HUMAN_VECTOR_INTERVENTIONS
#define OM_INTERVENTIONS_HUMAN_VECTOR_INTERVENTIONS

#include "interventions/Interventions.h"
#include <boost/ptr_container/ptr_list.hpp>

namespace OM { namespace interventions {
class HumanVectorInterventionParams;

/** A base class for interventions affecting human-vector interaction. */
class HumanVectorIntervention {
public:
    /** Deploy an intervention. */
    virtual void deploy( const HumanVectorInterventionParams& params ) =0;
    
    /// Get effect of deterrencies of interventions, as an attractiveness multiplier.
    virtual double relativeAttractiveness(const HumanInterventionEffect& params, size_t speciesIndex) const =0;
    /** Get the killing effect on mosquitoes before they've eaten as a survival
     * multiplier. */
    virtual double preprandialSurvivalFactor(const HumanInterventionEffect& params, size_t speciesIndex) const =0;
    /** Get the killing effect on mosquitoes after they've eaten as a survival
     * multiplier. */
    virtual double postprandialSurvivalFactor(const HumanInterventionEffect& params, size_t speciesIndex) const =0;
    
    /// Index of effect describing the intervention
    inline size_t getIndex() const { return index; }
    
    /// Checkpointing (write only)
    void operator& (ostream& stream) {
        index & stream;
        checkpoint( stream );
    }
    
protected:
    /// Set the effect index
    explicit HumanVectorIntervention( size_t index ) : index(index) {}
    
    /// Checkpointing: write
    virtual void checkpoint( ostream& stream ) =0;
private:
    size_t index;
};

/** A base class for human vector intervention parameters. */
class HumanVectorInterventionParams : public HumanInterventionEffect {
public:
    /** Create a new object to store human-specific details of deployment.
     * 
     * NOTE: no information about the target human is provided here; in theory
     * it can be provided if necessary. */
    virtual HumanVectorIntervention* makeHumanPart() const =0;
    virtual HumanVectorIntervention* makeHumanPart( istream& stream, size_t index ) const =0;
protected:
    explicit HumanVectorInterventionParams(size_t index) : HumanInterventionEffect(index) {}
};

/** Class to manage a set of vector interventions deployed to a human. */
class HumanVectorInterventions {
public:
    /** Deploy an intervention, provided with a certain parameters file. */
    void deploy( const HumanVectorInterventionParams& params );
    
    /** Get effect of deterrencies of interventions, as an attractiveness multiplier.
     * 
     * This is the product of (1 - deterrency) across all active interventions,
     * or 1 if no interventions are active.*/
    double relativeAttractiveness( size_t speciesIndex )const;
    /** Get the killing effect on mosquitoes before they've eaten as a survival
     * multiplier.
     * 
     * This is the product of (1 - pre_prandial_kill_factor) across all active
     * interventions, or 1 if no interventions are active. */
    double preprandialSurvivalFactor( size_t speciesIndex )const;
    /** Get the killing effect on mosquitoes after they've eaten as a survival
     * multiplier.
     * 
     * This is the product of (1 - post_prandial_kill_factor) across all active
     * interventions, or 1 if no interventions are active. */
    /// See ComponentParams::effect for a more detailed description.
    double postprandialSurvivalFactor( size_t speciesIndex )const;
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        checkpoint( stream );
    }
    
private:
    void checkpoint( ostream& stream );
    void checkpoint( istream& stream );
    
    typedef boost::ptr_list<HumanVectorIntervention> ActiveList;
    ActiveList active;
};

} }

#endif
