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

// Note: this file is included by exactly one source (InterventionManager.cpp)
// and contains definitions as well as declarations.

// The includes here are more for documentation than required.
#include "Global.h"
#include "util/errors.h"
#include "Monitoring/Survey.h"
#include "Clinical/ClinicalModel.h"
#include "Population.h"
#include "Transmission/TransmissionModel.h"
#include "util/random.h"
#include <schema/interventions.h>

namespace OM { namespace interventions {

// ———  TimedDeployment and derivatives  ———

/** Interface for timed deployment of an intervention. */
class TimedDeployment {
public:
    /// Create, passing time of deployment
    explicit TimedDeployment(SimTime deploymentTime) :
            time( deploymentTime )
    {
        if( deploymentTime < sim::zero() ){
            throw util::xml_scenario_error("timed intervention deployment: may not be negative");
        }else if( deploymentTime >= Monitoring::Survey::getLastSurveyTime() ){
            //TODO(date): output as a date
            cerr << "Warning: timed intervention deployment at time "<<(deploymentTime / sim::oneTS());
            cerr << " happens after last survey" << endl;
        }
    }
    virtual ~TimedDeployment() {}
    
    bool operator< (const TimedDeployment& that) const{
        return this->time < that.time;
    }
    
    virtual void deploy (OM::Population&) =0;
    
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const =0;
#endif
    
    // Read access required in this file; don't really need protection:
    SimTime time;
};

class DummyTimedDeployment : public TimedDeployment {
public:
    DummyTimedDeployment() :
        TimedDeployment( sim::zero() )
    {
        // TimedDeployment's ctor checks that the deployment time-step is
        // within the intervention period. We want this time to be after the
        // last time-step, so set the time here after TimedDeployment's ctor
        // check has been done (hacky).
        time = sim::future();
    }
    virtual void deploy (OM::Population&) {}
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const{
        out << time << "\t\t\t\t\tdummy (no interventions)";
    }
#endif
};

class TimedChangeHSDeployment : public TimedDeployment {
public:
    TimedChangeHSDeployment( const scnXml::ChangeHS::TimedDeploymentType& hs ) :
        TimedDeployment( sim::fromTS( hs.getTime() ) /*FIXME(schema): units*/ ),
        newHS( hs._clone() )
    {}
    virtual void deploy (OM::Population& population) {
        Clinical::ClinicalModel::changeHS( *newHS );
        delete newHS;
        newHS = 0;
    }
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const{
        out << time << "\t\t\t\t\tchange HS";
    }
#endif
    
private:
    scnXml::HealthSystem *newHS;
};

class TimedChangeEIRDeployment : public TimedDeployment {
public:
    TimedChangeEIRDeployment( const scnXml::ChangeEIR::TimedDeploymentType& nv ) :
        TimedDeployment(sim::fromTS( nv.getTime() ) /*FIXME(schema): units*/ ),
        newEIR( nv._clone() )
    {}
    virtual void deploy (OM::Population& population) {
        population.transmissionModel().changeEIRIntervention( *newEIR );
        delete newEIR;
        newEIR = 0;
    }
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const{
        out << time << "\t\t\t\t\tchange EIR";
    }
#endif
    
private:
    scnXml::NonVector *newEIR;
};

class TimedUninfectVectorsDeployment : public TimedDeployment {
public:
    TimedUninfectVectorsDeployment( SimTime deployTime ) :
        TimedDeployment( deployTime )
    {}
    virtual void deploy (OM::Population& population) {
        population.transmissionModel().uninfectVectors();
    }
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const{
        out << time << "\t\t\t\t\tuninfect vectors";
    }
#endif
};

#if 0
class TimedR_0Deployment : public TimedDeployment {
public:
    TimedR_0Deployment( SimTime deployTime ) :
        TimedDeployment( deployTime )
    {}
    virtual void deploy (OM::Population& population) {
        int i = (int)std::floor (util::random::uniform_01() * population.size());        // pick a human
        Population::Iter it = population.begin();
        while (i > 0){  // find human (can't use population[i])
            ++it;
            --i;
        }
        assert( i == 0 );
        assert( it != population.end() );
        it->R_0Vaccines();
        it->addInfection();
    }
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const{
        out << time << "\t\t\t\t\tR_0 special";
    }
#endif
};
#endif

void VaccineLimits::set( const scnXml::DeploymentBase& deploy ){
    if( deploy.getVaccMinPrevDoses().present() ){
        int v = deploy.getVaccMinPrevDoses().get();
        if( v < 0 ) throw util::xml_scenario_error( "vaccMinPrevDoses: min value is 0" );
        minPrevDoses = v;
    }
    if( deploy.getVaccMaxCumDoses().present() ){
        int v = deploy.getVaccMaxCumDoses().get();
        if( v < 0 ) throw util::xml_scenario_error( "vaccMaxCumDoses: min value is 0" );
        maxCumDoses = v;
    }
}

/// Base class for TimedHumanDeployment and ContinuousHumanDeployment
class HumanDeploymentBase {
protected:
    /**
     * @param deploy XML element describing deployment
     * @param intervention The intervention to deploy (list of components)
     * @param subPop Either ComponentId_pop or a sub-population to which deployment is restricted
     * @param complement Whether to take the complement of the sub-population
     *  to which deployment will be restricted
     */
    HumanDeploymentBase( const scnXml::DeploymentBase& deploy,
                         const HumanIntervention* intervention,
                         ComponentId subPop, bool complement ) :
            coverage( deploy.getCoverage() ),
            subPop( subPop ),
            complement( complement ),
            intervention( intervention )
    {
        if( !(coverage >= 0.0 && coverage <= 1.0) ){
            throw util::xml_scenario_error("intervention deployment coverage must be in range [0,1]");
        }
        vaccLimits.set( deploy );
    }
    
    inline void deployToHuman( Host::Human& human, Deployment::Method method ) const{
        intervention->deploy( human, method, vaccLimits );
    }
    
    double coverage;    // proportion coverage within group meeting above restrictions
    VaccineLimits vaccLimits;
    ComponentId subPop;      // ComponentId_pop if deployment is not restricted to a sub-population
    bool complement;
    const HumanIntervention *intervention;
};

/// Timed deployment of human-specific interventions
class TimedHumanDeployment : public TimedDeployment, protected HumanDeploymentBase {
public:
    /** 
     * @param mass XML element specifying the age range and compliance
     * (proportion of eligible individuals who receive the intervention).
     * @param intervention The HumanIntervention to deploy.
     * @param subPop Either ComponentId_pop or a sub-population to which deployment is restricted
     */
    TimedHumanDeployment( const scnXml::Mass& mass,
                           const HumanIntervention* intervention,
                           ComponentId subPop, bool complement ) :
        TimedDeployment( sim::fromTS( mass.getTime() ) /*FIXME(schema): units*/ ),
        HumanDeploymentBase( mass, intervention, subPop, complement ),
        minAge( sim::fromYearsN( mass.getMinAge() ) ),
        maxAge( sim::future() )
    {
        if( mass.getMaxAge().present() )
            maxAge = sim::fromYearsN( mass.getMaxAge().get() );
            
        if( minAge < sim::zero() || maxAge < minAge ){
            throw util::xml_scenario_error("timed intervention must have 0 <= minAge <= maxAge");
        }
    }
    
    virtual void deploy (OM::Population& population) {
        for (Population::Iter iter = population.begin(); iter != population.end(); ++iter) {
            SimTime age = iter->getAge();
            if( age >= minAge && age < maxAge ){
                if( subPop == interventions::ComponentId_pop || (iter->isInSubPop( subPop ) != complement) ){
                    if( util::random::bernoulli( coverage ) ){
                        deployToHuman( *iter, Deployment::TIMED );
                    }
                }
            }
        }
    }
    
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const{
        out << time << '\t'
            << minAge << '\t' << maxAge << '\t';
        if( subPop == ComponentId_pop ) out << "(none)";
        else out << subPop.id;
        out << '\t' << complement << '\t' << coverage << '\t';
        intervention->print_details( out );
    }
#endif
    
protected:
    // restrictions on deployment
    SimTime minAge, maxAge;
};

/// Timed deployment of human-specific interventions in cumulative mode
class TimedCumulativeHumanDeployment : public TimedHumanDeployment {
public:
    /** 
     * @param mass XML element specifying the age range and compliance
     * (proportion of eligible individuals who receive the intervention).
     * @param intervention The HumanIntervention to deploy.
     * @param subPop Either ComponentId_pop or a sub-population to which deployment is restricted
     * @param cumCuvId Id of component to test coverage for
     */
    TimedCumulativeHumanDeployment( const scnXml::Mass& mass,
                           const HumanIntervention* intervention,
                           ComponentId subPop, bool complement,
                           ComponentId cumCuvId ) :
        TimedHumanDeployment( mass, intervention, subPop, complement ),
        cumCovInd( cumCuvId )
    {
    }
    
    virtual void deploy (OM::Population& population) {
        // Cumulative case: bring target group's coverage up to target coverage
        vector<Host::Human*> unprotected;
        size_t total = 0;       // number of humans within age bound and optionally subPop
        for (Population::Iter iter = population.begin(); iter != population.end(); ++iter) {
            SimTime age = iter->getAge();
            if( age >= minAge && age < maxAge ){
                if( subPop == interventions::ComponentId_pop || (iter->isInSubPop( subPop ) != complement) ){
                    total+=1;
                    if( !iter->isInSubPop(cumCovInd) )
                        unprotected.push_back( &*iter );
                }
            }
        }
        
        if( total == 0 ) return;        // no humans to deploy to; avoid divide by zero
        double propProtected = static_cast<double>( total - unprotected.size() ) / static_cast<double>( total );
        if( propProtected < coverage ){
            // Proportion propProtected are already covered, so need to
            // additionally cover the proportion (coverage - propProtected),
            // selected from the list unprotected.
            double additionalCoverage = (coverage - propProtected) / (1.0 - propProtected);
            cerr << "cum deployment: prop protected " << propProtected << "; additionalCoverage " << additionalCoverage << "; total " << total << endl;
            for (vector<Host::Human*>::iterator iter = unprotected.begin();
                 iter != unprotected.end(); ++iter)
            {
                if( util::random::uniform_01() < additionalCoverage ){
                    deployToHuman( **iter, Deployment::TIMED );
                }
            }
        }
    }
    
protected:
    ComponentId cumCovInd;
};

class TimedVectorDeployment : public TimedDeployment {
public:
    TimedVectorDeployment( SimTime deployTime, size_t instance ) :
        TimedDeployment( deployTime ),
        inst(instance)
    {}
    virtual void deploy (OM::Population& population) {
      population.transmissionModel().deployVectorPopInterv(inst);
    }
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const{
        out << time << "\t\t\t\t\tvector";
    }
#endif
private:
    size_t inst;
};

// ———  ContinuousHumanDeployment  ———

/** Interface for continuous deployment of an intervention. */
class ContinuousHumanDeployment : protected HumanDeploymentBase {
public:
    /// Create, passing deployment age
    ContinuousHumanDeployment( const ::scnXml::ContinuousDeployment& elt,
                                 const HumanIntervention* intervention,
                                 ComponentId subPop, bool complement ) :
            HumanDeploymentBase( elt, intervention, subPop, complement ),
            //FIXME(schema): time units
            begin( sim::fromTS(elt.getBegin()) ),
            end( sim::fromTS(elt.getEnd()) ),
            deployAge( sim::fromYearsN( elt.getTargetAgeYrs() ) )
    {
        if( begin < sim::zero() || end < begin ){
            throw util::xml_scenario_error("continuous intervention must have 0 <= begin <= end");
        }
        if( deployAge <= sim::zero() ){
            ostringstream msg;
            msg << "continuous intervention with target age "<<elt.getTargetAgeYrs();
            msg << " years corresponds to time step "<<deployAge;
            msg << "; must be at least 1.";
            throw util::xml_scenario_error( msg.str() );
        }
        if( deployAge > sim::maxHumanAge() ){
            ostringstream msg;
            msg << "continuous intervention must have target age no greater than ";
            msg << sim::maxHumanAge().inYears();
            throw util::xml_scenario_error( msg.str() );
        }
    }
    
    /// For sorting
    inline bool operator<( const ContinuousHumanDeployment& that )const{
        return this->deployAge < that.deployAge;
    }
    
    /** Apply filters and potentially deploy.
     * 
     * @returns false iff this deployment (and thus all later ones in the
     *  ordered list) happens in the future. */
    bool filterAndDeploy( Host::Human& human, const Population& population ) const{
        SimTime age = human.getAge();
        if( deployAge > age ){
            // stop processing continuous deployments for this
            // human for now because remaining ones happen in the future
            return false;
        }else if( deployAge == age ){
            if( begin <= sim::intervNow() && sim::intervNow() < end &&
                ( subPop == interventions::ComponentId_pop || (human.isInSubPop( subPop ) != complement) ) &&
                util::random::uniform_01() < coverage )     // RNG call should be last test
            {
                deployToHuman( human, Deployment::CTS );
            }
        }//else: for some reason, a deployment age was missed; ignore it
        return true;
    }
    
#ifdef WITHOUT_BOINC
    inline void print_details( std::ostream& out )const{
        out << begin << '\t';
        if( end == sim::future() ) out << "(none)";
        else out << end;
        out << '\t' << deployAge << '\t';
        if( subPop == ComponentId_pop ) out << "(none)";
        else out << subPop.id;
        out << '\t' << complement << '\t' << coverage << '\t';
        intervention->print_details( out );
    }
#endif
    
protected:
    SimTime begin, end;    // first time step active and first time step no-longer active
    SimTime deployAge;
};

} }
