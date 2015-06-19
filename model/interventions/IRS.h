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

#ifndef OM_INTERVENTIONS_IRS
#define OM_INTERVENTIONS_IRS

#include "util/DecayFunction.h"
#include "Transmission/PerHost.h"
#include "schema/interventions.h"
#include <boost/shared_ptr.hpp>

namespace OM {
namespace interventions {
    using util::DecayFunction;
    using util::DecayFuncHet;
    using util::NormalSampler;
    using util::LognormalSampler;
    using Transmission::PerHostInterventionData;

class IRSComponent : public Transmission::HumanVectorInterventionComponent {
public:
    IRSComponent( ComponentId id, const scnXml::IRSDescription& elt,
        const map<string,size_t>& species_name_map );
    
    virtual void deploy( Host::Human& human, mon::Deploy::Method method, VaccineLimits )const;
    
    virtual Component::Type componentType() const;
    
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const;
#endif
    
    virtual PerHostInterventionData* makeHumanPart() const;
    virtual PerHostInterventionData* makeHumanPart( istream& stream, ComponentId ) const;
    
private:
    /** Per mosquito-species parameters for extended IRS model. */
    class IRSAnopheles {
    public:
        IRSAnopheles() :
            proportionProtected( numeric_limits<double>::signaling_NaN() ),
            proportionUnprotected( numeric_limits<double>::signaling_NaN() )
        {}
        void init(const scnXml::IRSDescription::AnophelesParamsType& elt,
                 double proportionUse, double maxInsecticide);
        
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
        
        /// Return x*proportionProtected + proportionUnprotected
        inline double byProtection(double x) const{
            return x*proportionProtected + proportionUnprotected;
        }
        
    private:
        class RelativeAttractiveness {
        public:
            RelativeAttractiveness();
            
            /** Set parameters.
            * 
            * It is checked that input parameters lie in a range such that
            * the relative availability is always in the range (0,1] — that is,
            * the deterrent can never be perfect, but can have zero effect. */
            void init(const scnXml::IRSDeterrency& elt);
            
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
            
            /** Set parameters.
            * 
            * It is checked that parameters lie in a suitible range, giving a
            * survival factor between 0 and 1. */
            void init(const scnXml::IRSKillingEffect& elt, bool postPrandial, double maxInsecticide);
            
            /** Calculate additional survival factor imposed by IRS on pre-/post-
            * prandial killing. Should be bounded to [0,1] and tend to 1 as the
            * IRS ages. */
            double survivalFactor( double insecticideContent )const;
        private:
            double BF,PF;   // base and insecticide factors
            double insecticideScaling;
            double invBaseSurvival; // stored for performance only
        };
        double proportionProtected;
        double proportionUnprotected;
        RelativeAttractiveness _relativeAttractiveness;
        SurvivalFactor _preprandialKillingEffect;
        SurvivalFactor _postprandialKillingEffect;
        
        friend class HumanIRS;
    };
    
    double sampleInitialInsecticide() const;
    
    NormalSampler initialInsecticide;
    double maxInsecticide;              // maximum initial insecticide
    boost::shared_ptr<DecayFunction> insecticideDecay;
    vector<IRSAnopheles> species; // vector specific params
    
    // This is sparse vector: only indexes corresponding to a IRS component are used
    // No memory management
    static vector<IRSComponent*> componentsByIndex;
    
    friend class HumanIRS;
};

/** Extended IRS model by OB and original model.
 * 
 * Each instance describes the effects of indoor residual spraying.
 */
class HumanIRS : public PerHostInterventionData {
public:
    HumanIRS( const IRSComponent& params );
    HumanIRS( istream& stream, ComponentId id );
    
    virtual void redeploy(const Transmission::HumanVectorInterventionComponent& params);
    
    /// Get remaining insecticide content based on initial amount and decay.
    inline double getInsecticideContent(const IRSComponent& params)const{
        SimTime age = sim::nowOrTs1() - deployTime;  // implies age 1 TS on first use
        double effectSurvival = params.insecticideDecay->eval( age,
                                              insecticideDecayHet );
        return initialInsecticide * effectSurvival;
    }
    
    /// Call once per time step to update holes
    virtual void update(Host::Human& human);
    
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
    // this is sampled for each deployment: initial insecticide content doesn't
    // depend on handling by the recipient
    double initialInsecticide;	// units: mg/m²
    
    // this parameter is sampled on first deployment, but never resampled for the same human:
    DecayFuncHet insecticideDecayHet;
};

} }

#endif
