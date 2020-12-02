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
#include "Clinical/ClinicalModel.h"
#include "Population.h"
#include "Transmission/TransmissionModel.h"
#include "util/random.h"
#include <schema/interventions.h>

namespace OM { namespace interventions {

struct ByDeployTime;    // forward decl for friend

// ———  TimedDeployment and derivatives  ———

/** Interface for timed deployment of an intervention. */
class TimedDeployment {
public:
    /// Create, passing time of deployment
    explicit TimedDeployment(SimDate deployDate) :
            date( deployDate )
    {
        if( deployDate < sim::startDate() ){
            throw util::xml_scenario_error("timed intervention deployment: may not be negative");
        }else if( deployDate >= sim::endDate() ){
            cerr << "Warning: timed intervention deployment at date "<< deployDate;
            cerr << " happens after last survey" << endl;
        }
    }
    virtual ~TimedDeployment() {}
    
    virtual void deploy (Population& population, Transmission::TransmissionModel& transmission) =0;
    
    virtual void print_details( std::ostream& out )const =0;
    
    // Read access required in this file; don't really need protection:
    SimDate date;
};

class DummyTimedDeployment : public TimedDeployment {
public:
    DummyTimedDeployment() :
        TimedDeployment( sim::startDate() /* hack */ )
    {
        // TimedDeployment's ctor checks that the deployment time-step is
        // within the intervention period. We want this time to be after the
        // last time-step, so set the time here after TimedDeployment's ctor
        // check has been done (hacky).
        date = SimDate::future();
    }
    virtual void deploy (Population& population, Transmission::TransmissionModel& transmission) {}
    virtual void print_details( std::ostream& out )const{
        out << date << "\t\t\t\t\tdummy (no interventions)";
    }
};

class TimedChangeHSDeployment : public TimedDeployment {
public:
    TimedChangeHSDeployment( SimDate date, const scnXml::ChangeHS::TimedDeploymentType& hs ) :
        TimedDeployment( date ),
        newHS( hs._clone() )
    {}
    virtual void deploy (Population& population, Transmission::TransmissionModel& transmission) {
        Clinical::ClinicalModel::setHS( *newHS );
        delete newHS;
        newHS = 0;
    }
    virtual void print_details( std::ostream& out )const{
        out << date << "\t\t\t\t\tchange HS";
    }
    
private:
    scnXml::HealthSystem *newHS;
};

class TimedChangeEIRDeployment : public TimedDeployment {
public:
    TimedChangeEIRDeployment( SimDate date, const scnXml::ChangeEIR::TimedDeploymentType& nv ) :
        TimedDeployment( date ),
        newEIR( nv._clone() )
    {}
    virtual void deploy (Population& population, Transmission::TransmissionModel& transmission) {
        transmission.changeEIRIntervention( *newEIR );
        delete newEIR;
        newEIR = 0;
    }
    virtual void print_details( std::ostream& out )const{
        out << date << "\t\t\t\t\tchange EIR";
    }
    
private:
    scnXml::NonVector *newEIR;
};

class TimedUninfectVectorsDeployment : public TimedDeployment {
public:
    TimedUninfectVectorsDeployment( SimDate date ) :
        TimedDeployment( date )
    {}
    virtual void deploy (Population& population, Transmission::TransmissionModel& transmission) {
        transmission.uninfectVectors();
    }
    virtual void print_details( std::ostream& out )const{
        out << date << "\t\t\t\t\tuninfect vectors";
    }
};

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
     * @param subPop Either ComponentId::wholePop() or a sub-population to which deployment is restricted
     * @param complement Whether to take the complement of the sub-population
     *  to which deployment will be restricted
     */
    HumanDeploymentBase( const scnXml::DeploymentBase& deploy,
                         shared_ptr<const HumanIntervention> intervention,
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
    
    inline void deployToHuman( Host::Human& human, mon::Deploy::Method method ) const{
        intervention->deploy( human, method, vaccLimits );
    }
    
    double coverage;    // proportion coverage within group meeting above restrictions
    VaccineLimits vaccLimits;
    ComponentId subPop;      // ComponentId::wholePop() if deployment is not restricted to a sub-population
    bool complement;
    shared_ptr<const HumanIntervention> intervention;
};

/// Timed deployment of human-specific interventions
class TimedHumanDeployment : public TimedDeployment, protected HumanDeploymentBase {
public:
    /** 
     * @param date Time of deployment
     * @param mass XML element specifying the age range and compliance
     * (proportion of eligible individuals who receive the intervention).
     * @param intervention The HumanIntervention to deploy.
     * @param subPop Either ComponentId::wholePop() or a sub-population to which deployment is restricted
     */
    TimedHumanDeployment( SimDate date,
                           const scnXml::MassDeployment& mass,
                           shared_ptr<const HumanIntervention> intervention,
                           ComponentId subPop, bool complement ) :
        TimedDeployment( date ),
        HumanDeploymentBase( mass, intervention, subPop, complement ),
        minAge( SimTime::fromYearsN( mass.getMinAge() ) ),
        maxAge( SimTime::future() )
    {
        if( mass.getMaxAge().present() )
            maxAge = SimTime::fromYearsN( mass.getMaxAge().get() );
            
        if( minAge < SimTime::zero() || maxAge < minAge ){
            throw util::xml_scenario_error("timed intervention must have 0 <= minAge <= maxAge");
        }
    }
    
    virtual void deploy (Population& population, Transmission::TransmissionModel& transmission) {
        for(Human& human : population) {
            SimTime age = human.age(sim::now());
            if( age >= minAge && age < maxAge ){
                if( subPop == ComponentId::wholePop() || (human.isInSubPop( subPop ) != complement) ){
                    if( human.rng().bernoulli( coverage ) ){
                        deployToHuman( human, mon::Deploy::TIMED );
                    }
                }
            }
        }
    }
    
    virtual void print_details( std::ostream& out )const{
        out << date << "\t"
            << minAge.inYears() << "y\t" << maxAge.inYears() << "t\t";
        if( subPop == ComponentId::wholePop() ) out << "(none)";
        else out << subPop.id;
        out << '\t' << complement << '\t' << coverage << '\t';
        intervention->print_details( out );
    }
    
protected:
    // restrictions on deployment
    SimTime minAge, maxAge;
};

/// Timed deployment of human-specific interventions in cumulative mode
class TimedCumulativeHumanDeployment : public TimedHumanDeployment {
public:
    /** 
     * @param date Time of deployment
     * @param mass XML element specifying the age range and compliance
     * (proportion of eligible individuals who receive the intervention).
     * @param intervention The HumanIntervention to deploy.
     * @param subPop Either ComponentId::wholePop() or a sub-population to which deployment is restricted
     * @param cumCuvId Id of component to test coverage for
     */
    TimedCumulativeHumanDeployment( SimDate date,
                           const scnXml::MassDeployment& mass,
                           shared_ptr<const HumanIntervention> intervention,
                           ComponentId subPop, bool complement,
                           ComponentId cumCuvId ) :
        TimedHumanDeployment( date, mass, intervention, subPop, complement ),
        cumCovInd( cumCuvId )
    {
    }
    
    virtual void deploy (Population& population, Transmission::TransmissionModel& transmission) {
        // Cumulative case: bring target group's coverage up to target coverage
        vector<Host::Human*> unprotected;
        size_t total = 0;       // number of humans within age bound and optionally subPop
        for(Population::Iter iter = population.begin(); iter != population.end(); ++iter) {
            SimTime age = iter->age(sim::now());
            if( age >= minAge && age < maxAge ){
                if( subPop == ComponentId::wholePop() || (iter->isInSubPop( subPop ) != complement) ){
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
            for(Human* human : unprotected) {
                if( human->rng().uniform_01() < additionalCoverage ){
                    deployToHuman( *human, mon::Deploy::TIMED );
                }
            }
        }
    }
    
protected:
    ComponentId cumCovInd;
};

class TimedVectorDeployment : public TimedDeployment {
public:
    TimedVectorDeployment( SimDate date, size_t instance ) :
        TimedDeployment( date ),
        inst(instance)
    {}
    virtual void deploy (Population& population, Transmission::TransmissionModel& transmission) {
      transmission.deployVectorPopInterv(inst);
    }
    virtual void print_details( std::ostream& out )const{
        out << date << "\t\t\t\t\tvector";
    }
    
private:
    size_t inst;
};

class TimedTrapDeployment : public TimedDeployment {
public:
    TimedTrapDeployment( SimDate date, size_t instance, double ratio, SimTime lifespan ) :
        TimedDeployment(date), inst(instance), ratio(ratio), lifespan(lifespan)
    {}
    virtual void deploy (Population& population, Transmission::TransmissionModel& transmission) {
        double number = population.size() * ratio;
        transmission.deployVectorTrap(inst, number, lifespan);
    }
    virtual void print_details( std::ostream& out )const{
        out << date << "\t\t\t\t\tvector trap";
    }
    
private:
    size_t inst;
    double ratio;
    SimTime lifespan;
};

class TimedNonHumanHostsDeployment : public TimedDeployment {
public:
    TimedNonHumanHostsDeployment( SimDate date, size_t instance, string name ) :
        TimedDeployment( date ),
        instance(instance),
        name(name)
    {}
    virtual void deploy (Population& population, Transmission::TransmissionModel& transmission) {
      transmission.deployNonHumanHostsInterv(instance, name);
    }
    virtual void print_details( std::ostream& out )const{
        out << date << "\t\t\t\t\tnhh";
    }
    
private:
    size_t instance;
    string name;
};

class TimedAddNonHumanHostsDeployment : public TimedDeployment {
public:
    TimedAddNonHumanHostsDeployment( SimDate date, const string &name, SimTime lifespan) :
        TimedDeployment( date ), name(name), lifespan(lifespan)
    {}
    virtual void deploy (Population& population, Transmission::TransmissionModel& transmission) {
        double popSize = population.size();
        transmission.deployAddNonHumanHosts(name, popSize, lifespan);
    }
    virtual void print_details( std::ostream& out )const{
        out << date << "\t\t\t\t\tnhh";
    }
    
private:
    string name;
    SimTime lifespan;
};

// ———  ContinuousHumanDeployment  ———

/** Interface for continuous deployment of an intervention. */
class ContinuousHumanDeployment : protected HumanDeploymentBase {
public:
    /// Create, passing deployment age
    ContinuousHumanDeployment( SimDate begin, SimDate end,
                                 const ::scnXml::ContinuousDeployment& elt,
                                 shared_ptr<const HumanIntervention> intervention,
                                 ComponentId subPop, bool complement ) :
            HumanDeploymentBase( elt, intervention, subPop, complement ),
            begin( begin ), end( end ),
            deployAge( SimTime::fromYearsN( elt.getTargetAgeYrs() ) )
    {
        if( begin < sim::startDate() || end < begin ){
            throw util::xml_scenario_error("continuous intervention must have startDate <= begin <= end");
        }
        if( deployAge <= SimTime::zero() ){
            ostringstream msg;
            msg << "continuous intervention with target age "<<elt.getTargetAgeYrs();
            msg << " years corresponds to time step " << deployAge.inSteps();
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
    
    /** Apply filters and potentially deploy.
     * 
     * @returns false iff this deployment (and thus all later ones in the
     *  ordered list) happens in the future. */
    bool filterAndDeploy( Host::Human& human, const Population& population ) const{
        SimTime age = human.age(sim::now());
        if( deployAge > age ){
            // stop processing continuous deployments for this
            // human for now because remaining ones happen in the future
            return false;
        }else if( deployAge == age ){
            auto now = sim::intervDate();
            if( begin <= now && now < end &&
                ( subPop == ComponentId::wholePop() ||
                    (human.isInSubPop( subPop ) != complement)
                ) &&
                human.rng().uniform_01() < coverage )     // RNG call should be last test
            {
                deployToHuman( human, mon::Deploy::CTS );
            }
        }//else: for some reason, a deployment age was missed; ignore it
        return true;
    }
    
    inline void print_details( std::ostream& out )const{
        out << begin << "\t";
        if( end == SimDate::future() ) out << "(none)";
        else out << end << 't';
        out << '\t' << deployAge.inYears() << "y\t";
        if( subPop == ComponentId::wholePop() ) out << "(none)";
        else out << subPop.id;
        out << '\t' << complement << '\t' << coverage << '\t';
        intervention->print_details( out );
    }
    
protected:
    SimDate begin, end;    // first time step active and first time step no-longer active
    SimTime deployAge;
    
    friend ByDeployTime;
};

struct ByDeployTime {
    bool operator() (const unique_ptr<TimedDeployment>& a, const unique_ptr<TimedDeployment>& b) const{
        return a->date < b->date;
    }
    
    bool operator()(const ContinuousHumanDeployment& a, const ContinuousHumanDeployment& b) const{
        return a.deployAge < b.deployAge;
    }
} byDeployTime;

} }
