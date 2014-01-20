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

// Note: this file is included by exactly one source (Interventions.cpp)
// and contains definitions as well as declarations.

// The includes here are more for documentation than required.
#include "util/errors.h"
#include "util/TimeStep.h"
#include "Monitoring/Surveys.h"
#include <schema/interventions.h>
#include "Clinical/CaseManagementCommon.h"
#include "Population.h"
#include "Transmission/TransmissionModel.h"

namespace OM { namespace interventions {

// ———  TimedDeployment and derivatives  ———

/** Interface for timed deployment of an intervention. */
class TimedDeployment {
public:
    /// Create, passing time of deployment
    explicit TimedDeployment(TimeStep deploymentTime) :
            time( deploymentTime )
    {
        if( deploymentTime < TimeStep(0) ){
            throw util::xml_scenario_error("timed intervention deployment: may not be negative");
        }else if( deploymentTime >= Monitoring::Surveys.getFinalTimestep() ){
            cerr << "Warning: timed intervention deployment at time "<<deploymentTime.asInt();
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
    TimeStep time;
};

class DummyTimedDeployment : public TimedDeployment {
public:
    DummyTimedDeployment() :
        TimedDeployment( TimeStep(0) )
    {
        // TimedDeployment's ctor checks that the deployment time-step is
        // within the intervention period. We want this time to be after the
        // last time-step, so set the time here after TimedDeployment's ctor
        // check has been done (hacky).
        time = TimeStep::future;
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
        TimedDeployment( TimeStep( hs.getTime() ) ),
        newHS( hs._clone() )
    {}
    virtual void deploy (OM::Population& population) {
        Clinical::CaseManagementCommon::changeHealthSystem( *newHS );
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
        TimedDeployment( TimeStep( nv.getTime() ) ),
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
    TimedUninfectVectorsDeployment( TimeStep deployTime ) :
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

class TimedR_0Deployment : public TimedDeployment {
public:
    TimedR_0Deployment( TimeStep deployTime ) :
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

/// Timed deployment of human-specific interventions
class TimedHumanDeployment : public TimedDeployment {
public:
    /** 
     * @param mass XML element specifying the age range and compliance
     * (proportion of eligible individuals who receive the intervention).
     * @param intervention The HumanIntervention to deploy.
     * @param cohort The cohort to which to deploy, or max value */
    TimedHumanDeployment( const scnXml::Mass& mass,
                           const HumanIntervention* intervention,
                           size_t cohort ) :
        TimedDeployment( TimeStep( mass.getTime() ) ),
        minAge( TimeStep::fromYears( mass.getMinAge() ) ),
        maxAge( TimeStep::fromYears( mass.getMaxAge() ) ),
        cohort( cohort ),
        coverage( mass.getCoverage() ),
        intervention( intervention )
    {
        if( !(coverage >= 0.0 && coverage <= 1.0) ){
            throw util::xml_scenario_error("timed intervention coverage must be in range [0,1]");
        }
        if( minAge < TimeStep(0) || maxAge < minAge ){
            throw util::xml_scenario_error("timed intervention must have 0 <= minAge <= maxAge");
        }
    }
    
    virtual void deploy (OM::Population& population) {
        for (Population::Iter iter = population.begin(); iter != population.end(); ++iter) {
            TimeStep age = TimeStep::simulation - iter->getDateOfBirth();
            if( age >= minAge && age < maxAge ){
                if( iter->isInCohort( cohort) ){
                    if( util::random::uniform_01() < coverage ){
                        intervention->deploy( *iter, Deployment::TIMED );
                    }
                }
            }
        }
    }
    
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const{
        out << time << '\t'
            << minAge << '\t' << maxAge << '\t';
        if( cohort == numeric_limits<size_t>::max() ) out << "(none)";
        else out << cohort;
        out << '\t' << coverage << '\t';
        intervention->print_details( out );
    }
#endif
    
protected:
    // restrictions on deployment
    TimeStep minAge;
    TimeStep maxAge;
    size_t cohort;
    double coverage;    // proportion coverage within group meeting above restrictions
    const HumanIntervention *intervention;
};

/// Timed deployment of human-specific interventions in cumulative mode
class TimedCumulativeHumanDeployment : public TimedHumanDeployment {
public:
    /** 
     * @param mass XML element specifying the age range and compliance
     * (proportion of eligible individuals who receive the intervention).
     * @param intervention The HumanIntervention to deploy.
     * @param effect_index Index of effect to test coverage for
     * @param maxAge Maximum time-span to consider a deployed effect still to be effective */
    TimedCumulativeHumanDeployment( const scnXml::Mass& mass,
                           const HumanIntervention* intervention,
                           size_t cohort,
                           size_t effect_index, TimeStep maxAge ) :
        TimedHumanDeployment( mass, intervention, cohort ),
        cumCovInd( effect_index ), maxInterventionAge( maxAge )
    {
    }
    
    virtual void deploy (OM::Population& population) {
        // Cumulative case: bring target group's coverage up to target coverage
        vector<Host::Human*> unprotected;
        size_t total = 0;       // number of humans within age bound and optionally cohort
        for (Population::Iter iter = population.begin(); iter != population.end(); ++iter) {
            TimeStep age = TimeStep::simulation - iter->getDateOfBirth();
            if( age >= minAge && age < maxAge ){
                if( iter->isInCohort( cohort ) ){
                    total+=1;
                    if( iter->needsRedeployment(cumCovInd, maxInterventionAge) )
                        unprotected.push_back( &*iter );
                }
            }
        }
        
        double propProtected = static_cast<double>( total - unprotected.size() ) / static_cast<double>( total );
        if( propProtected < coverage ){
            // Proportion propProtected are already covered, so need to
            // additionally cover the proportion (coverage - propProtected),
            // selected from the list unprotected.
            double additionalCoverage = (coverage - propProtected) / (1.0 - propProtected);
            for (vector<Host::Human*>::iterator iter = unprotected.begin();
                 iter != unprotected.end(); ++iter)
            {
                if( util::random::uniform_01() < additionalCoverage ){
                    intervention->deploy( **iter, Deployment::TIMED );
                }
            }
        }
    }
    
protected:
    size_t cumCovInd;
    // max age at which an intervention is considered not to need replacement
    TimeStep maxInterventionAge;
};

class TimedVectorDeployment : public TimedDeployment {
public:
    TimedVectorDeployment( TimeStep deployTime, size_t instance ) :
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

} }
