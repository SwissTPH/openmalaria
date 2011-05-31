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

#ifndef TRANSMISSION_ITN
#define TRANSMISSION_ITN

#include "util/DecayFunction.h"
#include "util/sampler.h"
#include "schema/interventions.h"

namespace OM { namespace Transmission {
    using util::DecayFunction;
    using util::DecayFuncHet;
    using util::NormalSampler;
    using util::LognormalSampler;

/** Constant parameters for extended ITN model. */
class ITNParams {
public:
    ITNParams() : ripFactor( numeric_limits<double>::signaling_NaN() ) {}
    /** Set parameters from elt. */
    double init( const scnXml::ITNDescription& elt);
    
private:
    NormalSampler initialInsecticide;
    LognormalSampler holeRate;
    LognormalSampler ripRate;
    double ripFactor;   // factor expressing how significant rips are in comparison to holes
    shared_ptr<DecayFunction> insecticideDecay;
    
    friend class ITN;
};

/** Per mosquito-species parameters for extended ITN model. */
class ITNAnophelesParams {
public:
    ITNAnophelesParams( const ITNParams* b ) :
        base( b ),
        proportionProtected( numeric_limits<double>::signaling_NaN() ),
        proportionUnprotected( numeric_limits<double>::signaling_NaN() )
    {}
    void init(const scnXml::ITNDescription::AnophelesParamsType& elt, double proportionUse);
    
    /// Get deterrency. See ComponentParams::effect for a more detailed description.
    inline double relativeAvailability( double holeIndex, double insecticideContent )const{
        return proportionProtected * _relativeAvailability.relativeAvailability( holeIndex, insecticideContent ) + proportionUnprotected;
    }
    /// Get killing effect on mosquitoes before feeding.
    /// See ComponentParams::effect for a more detailed description.
    inline double preprandialSurvivalFactor( double holeIndex, double insecticideContent )const{
        return proportionProtected * _preprandialKillingEffect.survivalFactor( holeIndex, insecticideContent ) + proportionUnprotected;
    }
    /// Get killing effect on mosquitoes after they've eaten.
    /// See ComponentParams::effect for a more detailed description.
    inline double postprandialSurvivalFactor( double holeIndex, double insecticideContent )const{
        return proportionProtected * _postprandialKillingEffect.survivalFactor( holeIndex, insecticideContent ) + proportionUnprotected;
    }
    
private:
    class RelativeAvailability {
    public:
        RelativeAvailability();
        
        /** Set parameters.
         * 
         * It is checked that input parameters lie in a range such that
         * the relative availability is always in the range (0,1] — that is,
         * the deterrent can never be perfect, but can have zero effect. */
        void init(const scnXml::ITNAvailEffect& elt);
        
        /** Calculate effect. Positive is interpreted as having a positive effect
        * (thus decreasing availability or survival) and negative as having a
        * negative effect. Effect is not bounded, though it tends to
        * zero as holeIndex becomes large and insecticideContent tends to zero,
        * and parameters should be defined such that it is always in the
        * range [0,1]. */
        double relativeAvailability( double holeIndex, double insecticideContent )const;
        
    protected:
        double lHF, lPF, lIF;      // logs of hole, insecticide and interaction factors
        double holeScaling, insecticideScaling;
    };
    class SurvivalFactor {
    public:
        SurvivalFactor();
        
        /** Set parameters.
         * 
         * It is checked that parameters lie in a suitible range, giving a
         * survival factor between 0 and 1. */
        void init(const scnXml::ITNKillingEffect& elt);
        
        /** Calculate additional survival factor imposed by nets on pre-/post-
         * prandial killing. Should be bounded to [0,1] and tend to 1 as the
         * net ages. */
        double survivalFactor( double holeIndex, double insecticideContent )const;
    private:
        double BF, HF, PF, IF;	// base, hole, insecticide and interaction factors
        double holeScaling, insecticideScaling;
        double invBaseSurvival; // stored for performance only
    };
    const ITNParams* base;
    double proportionProtected;
    double proportionUnprotected;
    RelativeAvailability _relativeAvailability;
    SurvivalFactor _preprandialKillingEffect;
    SurvivalFactor _postprandialKillingEffect;
    
    friend class ITN;
};

/** Extended ITN model by OB.
 * 
 * Each instance describes a hypothetical net (or no net).
 */
class ITN {
public:
    ITN () :
        holeIndex( numeric_limits<double>::signaling_NaN() ),
        initialInsecticide( numeric_limits<double>::signaling_NaN() ),
        holeRate( numeric_limits<double>::signaling_NaN() ),
        ripRate( numeric_limits<double>::signaling_NaN() )
    {
        nHoles = 0;
    }
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        deployTime & stream;
        nHoles & stream;
        holeIndex & stream;
        initialInsecticide & stream;
        holeRate & stream;
        ripRate & stream;
        insecticideDecayHet & stream;
    }
    
    void deploy(const ITNParams& params);
    inline TimeStep timeOfDeployment()const{
        return deployTime;
    }
    
    /// Call once per timestep to update holes
    void update(const ITNParams& params);
    
    /// Get deterrency. See ComponentParams::effect for a more detailed description.
    double relativeAvailability(const ITNAnophelesParams& params) const;
    /// Get killing effect on mosquitoes before they've eaten.
    /// See ComponentParams::effect for a more detailed description.
    double preprandialSurvivalFactor(const ITNAnophelesParams& params) const;
    /// Get killing effect on mosquitoes after they've eaten.
    /// See ComponentParams::effect for a more detailed description.
    double postprandialSurvivalFactor(const ITNAnophelesParams& params) const;
    
private:
    // these parameters express the current state of the net:
    TimeStep deployTime;	// time of deployment or TimeStep::never
    int nHoles;			// total number of holes
    double holeIndex;	// a measure of both the number and size of holes
    double initialInsecticide;	// TODO: units; mg/m²?
    
    // these parameters are sampled from log-normal per net, but thereafter constant:
    double holeRate;	// rate at which new holes are created
    double ripRate;		// rate at which holes are enlarged
    DecayFuncHet insecticideDecayHet;
};

} }

#endif
