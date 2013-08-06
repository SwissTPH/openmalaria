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

#ifndef OM_INTERVENTIONS_GVI
#define OM_INTERVENTIONS_GVI

#include "util/DecayFunction.h"
#include "util/sampler.h"
#include "schema/interventions.h"
#include <boost/shared_ptr.hpp>

namespace OM {
namespace Transmission {
    // forward declare:
    class TransmissionModel;
}
namespace interventions {
    using util::DecayFunction;
    using util::DecayFuncHet;
    using util::NormalSampler;
    using util::LognormalSampler;
    using boost::shared_ptr;

/** Constant parameters for generic vector intervention model. */
class GVIParams {
public:
    GVIParams() :  maxInsecticide(numeric_limits< double >::signaling_NaN()) {}
    /** Set parameters */
    void init( const scnXml::GVIDescription& elt);
    
private:
    NormalSampler initialInsecticide;
    double maxInsecticide;		// maximum initial insecticide
    shared_ptr<DecayFunction> decay;
    
    friend class GVI;
    friend class GVIAnophelesParams;
};

/** Per mosquito-species parameters for generic vector intervention model. */
class GVIAnophelesParams {
public:
    GVIAnophelesParams( const GVIParams* b ) :
        base( b ),
        proportionProtected( numeric_limits<double>::signaling_NaN() ),
        proportionUnprotected( numeric_limits<double>::signaling_NaN() ),
        _relativeAttractiveness( numeric_limits<double>::signaling_NaN() ),
        _preprandialKillingEffect( numeric_limits<double>::signaling_NaN() ),
        _postprandialKillingEffect( numeric_limits<double>::signaling_NaN() )
    {}
    void init(const GVIParams& params,
              const scnXml::GVIDescription::AnophelesParamsType& elt);
    
    /// Return x*proportionProtected + proportionUnprotected
    inline double byProtection(double x) const{
        return x*proportionProtected + proportionUnprotected;
    }
    
private:
    const GVIParams* base;
    double proportionProtected;
    double proportionUnprotected;
    double _relativeAttractiveness;
    double _preprandialKillingEffect;
    double _postprandialKillingEffect;
    
    friend class GVI;
};

/** Low-level (generic) vector intervention model. Has three effects:
 * deterrency, pre-prandial killing and post-prandial killing.
 * 
 * This is the per-host (but not per vector) part.
 */
class GVI {
public:
    GVI (const Transmission::TransmissionModel& tm);
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        deployTime & stream;
        initialInsecticide & stream;
        decayHet & stream;
    }
    
    void deploy(const GVIParams& params);
    inline TimeStep timeOfDeployment()const{
        return deployTime;
    }
    /** This is the survival factor of the effect. */
    inline double getEffectSurvival(const GVIParams& params)const{
        return params.decay->eval (TimeStep::simulation - deployTime,
                                              decayHet);
    }
    
    /// Get deterrency. See ComponentParams::effect for a more detailed description.
    double relativeAttractiveness(const GVIAnophelesParams& params) const;
    /// Get killing effect on mosquitoes before they've eaten.
    /// See ComponentParams::effect for a more detailed description.
    double preprandialSurvivalFactor(const GVIAnophelesParams& params) const;
    /// Get killing effect on mosquitoes after they've eaten.
    /// See ComponentParams::effect for a more detailed description.
    double postprandialSurvivalFactor(const GVIAnophelesParams& params) const;
    
private:
    // these parameters express the current state of the intervention:
    TimeStep deployTime;	// time of deployment or TimeStep::never
    double initialInsecticide;	// units: mg/m²
    
    // this parameters are sampled from log-normal per intervention, but thereafter constant:
    DecayFuncHet decayHet;
};

} }

#endif
