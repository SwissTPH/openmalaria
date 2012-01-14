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

#ifndef TRANSMISSION_IRS
#define TRANSMISSION_IRS

#include "util/DecayFunction.h"
#include "util/sampler.h"
#include "schema/interventions.h"
#include <boost/shared_ptr.hpp>

namespace OM { namespace Transmission {
    using util::DecayFunction;
    using util::DecayFuncHet;
    using util::NormalSampler;
    using util::LognormalSampler;
    using boost::shared_ptr;
    // forward declare:
    class TransmissionModel;

/** Constant parameters for extended IRS model. */
class IRSParams {
public:
    IRSParams() : simpleModel(false), maxInsecticide(numeric_limits< double >::signaling_NaN()) {}
    /** Set parameters for the new model from elt. */
    void init( const scnXml::IRSDescription& elt);
    /** Set parameters for the old model from elt. Don't call both! */
    void init( const scnXml::IRSSimpleDescription& elt);
    
private:
    // If true, use the older model with direct decay of effect; otherwise,
    // Use the Briet model with decay of insecticide (similar to ITN model).
    bool simpleModel;
    NormalSampler initialInsecticide;
    double maxInsecticide;		// maximum initial insecticide
    shared_ptr<DecayFunction> insecticideDecay;
    
    friend class IRS;
    friend class IRSAnophelesParams;
};

/** Per mosquito-species parameters for extended IRS model. */
class IRSAnophelesParams {
public:
    IRSAnophelesParams( const IRSParams* b ) :
        base( b )
    {}
    void init(const IRSParams& params,
              const scnXml::IRSDescription::AnophelesParamsType& elt);
    void init(const IRSParams& params,
              const scnXml::IRSSimpleDescription::AnophelesParamsType& elt);
    
    /// Get deterrency. See ComponentParams::effect for a more detailed description.
    inline double relativeAttractiveness( double insecticideContent )const{
        return _relativeAttractiveness.relativeAttractiveness( insecticideContent );
    }
    /// Get killing effect on mosquitoes before feeding.
    /// See ComponentParams::effect for a more detailed description.
    inline double preprandialSurvivalFactor( double insecticideContent )const{
        return _preprandialKillingEffect.survivalFactor( insecticideContent );
    }
    /// Get killing effect on mosquitoes after they've eaten.
    /// See ComponentParams::effect for a more detailed description.
    inline double postprandialSurvivalFactor( double insecticideContent )const{
        return _postprandialKillingEffect.survivalFactor( insecticideContent );
    }
    
private:
    class RelativeAttractiveness {
    public:
        RelativeAttractiveness();
        
        // for the old model: use lPF instead as the deterrency
        inline double oldDeterrency() const { return lPF; }
        inline void oldDeterrency(double d) { lPF = d; }
        
        /** Set parameters.
         * 
         * It is checked that input parameters lie in a range such that
         * the relative availability is always in the range (0,1] — that is,
         * the deterrent can never be perfect, but can have zero effect. */
        void init(const OM::Transmission::IRSParams& params,
                  const scnXml::IRSDeterrency& elt);
        
        /** Calculate effect. Positive is interpreted as having a positive effect
        * (thus decreasing availability or survival) and negative as having a
        * negative effect. Effect is not bounded, though it tends to
        * zero as holeIndex becomes large and insecticideContent tends to zero,
        * and parameters should be defined such that it is always in the
        * range [0,1]. */
        double relativeAttractiveness( double insecticideContent )const;
        
    protected:
        double lPF;      // log of insecticide factor
        double insecticideScaling;
    };
    class SurvivalFactor {
    public:
        SurvivalFactor();
        
        // for the old model: use PF instead as the effect
        inline double oldEffect() const { return PF; }
        inline void oldEffect(double e) { PF = e; }
        
        /** Set parameters.
         * 
         * It is checked that parameters lie in a suitible range, giving a
         * survival factor between 0 and 1. */
        void init(const OM::Transmission::IRSParams& params,
                  const scnXml::IRSKillingEffect& elt, bool postPrandial);
        
        /** Calculate additional survival factor imposed by IRS on pre-/post-
         * prandial killing. Should be bounded to [0,1] and tend to 1 as the
         * IRS ages. */
        double survivalFactor( double insecticideContent )const;
    private:
        double BF,PF;	// base and insecticide factors
        double insecticideScaling;
        double invBaseSurvival; // stored for performance only
    };
    const IRSParams* base;
    RelativeAttractiveness _relativeAttractiveness;
    SurvivalFactor _preprandialKillingEffect;
    SurvivalFactor _postprandialKillingEffect;
    
    friend class IRS;
};

/** Extended IRS model by OB and original model.
 * 
 * Each instance describes the effects of indoor residual spraying.
 */
class IRS {
public:
    IRS (const TransmissionModel& tm);
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        deployTime & stream;
        initialInsecticide & stream;
        insecticideDecayHet & stream;
    }
    
    void deploy(const IRSParams& params);
    inline TimeStep timeOfDeployment()const{
        return deployTime;
    }
    /** Old model: this is the survival factor of the effect. New model: not
     * used, except as part of getInsecticideContent below. */
    inline double getEffectSurvival(const IRSParams& params)const{
        return params.insecticideDecay->eval (TimeStep::simulation - deployTime,
                                              insecticideDecayHet);
    }
    /// Get remaining insecticide content based on initial amount and decay.
    inline double getInsecticideContent(const IRSParams& params)const{
        return initialInsecticide * getEffectSurvival( params );
    }
    
    /// Get deterrency. See ComponentParams::effect for a more detailed description.
    double relativeAttractiveness(const IRSAnophelesParams& params) const;
    /// Get killing effect on mosquitoes before they've eaten.
    /// See ComponentParams::effect for a more detailed description.
    double preprandialSurvivalFactor(const IRSAnophelesParams& params) const;
    /// Get killing effect on mosquitoes after they've eaten.
    /// See ComponentParams::effect for a more detailed description.
    double postprandialSurvivalFactor(const IRSAnophelesParams& params) const;
    
private:
    // these parameters express the current state of the IRS:
    TimeStep deployTime;	// time of deployment or TimeStep::never
    double initialInsecticide;	// units: mg/m²
    
    // these parameters are sampled from log-normal per IRS, but thereafter constant:
    //Old model: used as heterogeneity of general decay
    DecayFuncHet insecticideDecayHet;
};

} }

#endif
