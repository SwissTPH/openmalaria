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

#ifndef OM_INTERVENTIONS_GVI
#define OM_INTERVENTIONS_GVI

#include "Transmission/PerHost.h"
#include "interventions/Interfaces.hpp"
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
    GVIEffect( EffectId, const scnXml::GVIDescription& elt,
               const map<string,size_t>& species_name_map );
    
    void deploy( Host::Human& human, Deployment::Method method, VaccineLimits )const;
    
    virtual Effect::Type effectType() const;
    
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const;
#endif
    
    virtual PerHostInterventionData* makeHumanPart() const;
    virtual PerHostInterventionData* makeHumanPart( istream& stream, EffectId ) const;
    
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
    HumanGVI( istream& stream, EffectId );
    
    virtual void redeploy( const Transmission::HumanVectorInterventionEffect& );
    
    /** This is the survival factor of the effect. */
    inline double getEffectSurvival(const GVIEffect& params)const{
        return params.decay->eval (TimeStep::simulation - deployTime, decayHet);
    }
    
    virtual void update();
    
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
    // this parameter is sampled on first deployment, but never resampled for the same human:
    DecayFuncHet decayHet;
};

} }

#endif
