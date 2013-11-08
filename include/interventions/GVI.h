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

#include "Interventions.h"
#include "Transmission/PerHost.h"
#include "util/DecayFunction.h"
#include "util/sampler.h"
#include "schema/interventions.h"
#include <boost/shared_ptr.hpp>

namespace OM { namespace interventions {
    using util::DecayFunction;
    using util::DecayFuncHet;
    using util::NormalSampler;
    using util::LognormalSampler;
    using boost::shared_ptr;
    using Transmission::PerHostInterventionData;

/** Constant parameters for generic vector intervention model. */
class GVIEffect : public Transmission::HumanVectorInterventionEffect {
public:
    /** Initialise parameters
     * 
     * @param elt Effect description from XML
     * @param species_name_map Map of species names to indices.
     */
    GVIEffect( size_t index, const scnXml::GVIDescription& elt,
               const map<string,size_t>& species_name_map );
    
    void deploy( Host::Human& human, Deployment::Method method )const;
    
    virtual Effect::Type effectType() const;
    
    virtual PerHostInterventionData* makeHumanPart() const;
    virtual PerHostInterventionData* makeHumanPart( istream& stream, size_t index ) const;
    
private:
    /** Per mosquito-species parameters for generic vector intervention model. */
    class GVIAnopheles {
    public:
        GVIAnopheles() :
            proportionProtected( numeric_limits<double>::signaling_NaN() ),
            proportionUnprotected( numeric_limits<double>::signaling_NaN() ),
            deterrency( numeric_limits<double>::signaling_NaN() ),
            preprandialKilling( numeric_limits<double>::signaling_NaN() ),
            postprandialKilling( numeric_limits<double>::signaling_NaN() )
        {}
        void init(const scnXml::GVIDescription::AnophelesParamsType& elt);
        
        /// Return x*proportionProtected + proportionUnprotected
        inline double byProtection(double x) const{
            return x*proportionProtected + proportionUnprotected;
        }
        
        double proportionProtected;
        double proportionUnprotected;
        double deterrency;
        double preprandialKilling;
        double postprandialKilling;
        
        friend class HumanGVI;
    };
    
    NormalSampler initialInsecticide;
    shared_ptr<DecayFunction> decay;
    vector<GVIAnopheles> species;  // vector specific params
    
    // This is sparse vector: only indexes corresponding to a GVI effect are used
    // No memory management
    static vector<GVIEffect*> effectsByIndex;
    
    friend class HumanGVI;
};

/** Low-level (generic) vector intervention model. Has three effects:
 * deterrency, pre-prandial killing and post-prandial killing.
 * 
 * This is the per-host (but not per vector) part.
 */
class HumanGVI : public PerHostInterventionData {
public:
    HumanGVI( const GVIEffect& params );
    HumanGVI( istream& stream, size_t index );
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        deployTime & stream;
        initialInsecticide & stream;
        decayHet & stream;
    }
    
    virtual void redeploy();
    
    inline TimeStep timeOfDeployment()const{
        return deployTime;
    }
    /** This is the survival factor of the effect. */
    inline double getEffectSurvival(const GVIEffect& params)const{
        return params.decay->eval (TimeStep::simulation - deployTime, decayHet);
    }
    
    /// Get deterrency. See ComponentParams::effect for a more detailed description.
    virtual double relativeAttractiveness(size_t speciesIndex) const;
    /// Get killing effect on mosquitoes before they've eaten.
    /// See ComponentParams::effect for a more detailed description.
    virtual double preprandialSurvivalFactor(size_t speciesIndex) const;
    /// Get killing effect on mosquitoes after they've eaten.
    /// See ComponentParams::effect for a more detailed description.
    virtual double postprandialSurvivalFactor(size_t speciesIndex) const;
    
protected:
    virtual void checkpoint( ostream& stream );
    
private:
    // these parameters express the current state of the intervention:
    TimeStep deployTime;	// time of deployment or TimeStep::never
    double initialInsecticide;	// units: mg/m²
    
    // this parameters are sampled from log-normal per intervention, but thereafter constant:
    DecayFuncHet decayHet;
};

} }

#endif
