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

#include "Interventions.h"
#include "Host/Vaccine.h"
#include "Population.h"
#include "util/random.h"
#include "Clinical/ESCaseManagement.h"
#include "WithinHost/DescriptiveIPTWithinHost.h"

namespace OM {

// -----  TimedIntervention and derivatives  -----

class DummyIntervention : public TimedIntervention {
public:
    DummyIntervention() :
        TimedIntervention( TimeStep::future )
    {}
    virtual void deploy () {}
};

TimedIntervention::TimedIntervention(TimeStep deploymentTime) :
    time( deploymentTime )
{}

/// Deployment of mass-to-human interventions
class TimedMassIntervention : public TimedIntervention {
public:
    /** 
     * @param mass XML element specifying the age range and compliance
     * (proportion of eligible individuals who receive the intervention).
     * @param pop A reference to the population
     * @param intervention A member-function pointer to a "void func ()" function
     * within human which activates the intervention. */
    TimedMassIntervention( const scnXml::Mass& mass,
                           Population& pop,
                           void (Host::Human::*deployIntervention)() ) :
        TimedIntervention( TimeStep( mass.getTime() ) ),
        minAge( TimeStep::fromYears( mass.getMinAge() ) ),
        maxAge( TimeStep::fromYears( mass.getMaxAge() ) ),
        cohortOnly( mass.getCohort() ),
        coverage( mass.getCoverage() ),
        population( pop ),
        intervention( deployIntervention )
    {}
    
    virtual void deploy () {
        Population::HumanPop& popList = population.getList();
        for (Population::HumanIter iter = popList.begin(); iter != popList.end(); ++iter) {
            //TODO: note that this won't be exactly the same as before, so we may want to revert to old behaviour to check results exactly
            TimeStep age = TimeStep::simulation - iter->getDateOfBirth();
            if( age >= minAge && age < maxAge ){
                if( !cohortOnly || iter->getInCohort() ){
                    if( util::random::uniform_01() < coverage ){
                        // This is UGLY syntax. It just means call intervention() on the human pointed by iter.
                        ( (*iter).*intervention) ();
                    }
                }
            }
        }
    }
    
protected:
    // restrictions on deployment
    TimeStep minAge;
    TimeStep maxAge;
    bool cohortOnly;
    double coverage;    // proportion coverage within group meeting above restrictions
    Population& population;
    void (Host::Human::*intervention) ();       // callback: per-human deployment
};

/// Deployment of mass-to-human interventions with cumulative-deployment support
class TimedMassCumIntervention : public TimedMassIntervention {
public:
    /** As massIntervention, but supports "increase to target coverage" mode:
     * Deployment is only to unprotected humans and brings
     * total coverage up to the level given in description.
     * 
     * @param mass XML element specifying the age range and compliance
     * (proportion of eligible individuals who receive the intervention).
     * @param pop A reference to the population
     * @param intervention A member-function pointer to a "void func ()" function
     * within human which activates the intervention.
     * @param isProtectedCb A member-function pointer to a
     * "bool func (TimeStep maxAge)" function on a Human which returns true if
     * the Human is still protected by an intervention of the type in question
     * which is no older than maxAge. */
    TimedMassCumIntervention( const scnXml::MassCum& mass,
                              Population& pop,
                              void (Host::Human::*deployIntervention)(),
                              bool (Host::Human::*isProtectedCb) (TimeStep) const ) :
        TimedMassIntervention( mass, pop, deployIntervention ),
        isProtected( isProtectedCb ),
        maxInterventionAge( TimeStep::fromYears( mass.getCumulativeWithMaxAge().get() ) )
    {}
    
    void deploy(){
        /*FIXME
        if( mass.getCumulativeWithMaxAge().present() == false ){
            // Usual case: simply deploy to coverage% of target group
            massIntervention( mass, intervention );
            return;
        }
        */
        // Cumulative case: bring target group's coverage up to target coverage
        
        Population::HumanPop& popList = population.getList();
        vector<Host::Human*> unprotected;
        size_t total = 0;       // number of humans within age bound and optionally cohort
        for (Population::HumanIter iter = popList.begin(); iter != popList.end(); ++iter) {
            TimeStep age = TimeStep::simulation - iter->getDateOfBirth();
            if( age >= minAge && age < maxAge ){
                if( !cohortOnly || iter->getInCohort() ){
                    total+=1;
                    if( !((*iter).*isProtected)(maxInterventionAge) )
                        unprotected.push_back( &*iter );
                }
            }
        }
        
        double propProtected = static_cast<double>( total - unprotected.size() ) / static_cast<double>( total );
        if( propProtected < coverage ){
            // In range (0,1]:
            double additionalCoverage = (coverage - propProtected) / (1.0 - propProtected);
            for (Population::HumanIter iter = popList.begin(); iter != popList.end(); ++iter) {
                if( util::random::uniform_01() < additionalCoverage ){
                    ( (*iter).*intervention) ();
                }
            }
        }
    }
    
private:
    // callback to ascertain whether a human is still under protection from an
    // intervention young enough not to need replacement
    bool (Host::Human::*isProtected) (TimeStep) const;
    // max age at which an intervention is considered not to need replacement
    TimeStep maxInterventionAge;
};


// -----  InterventionManager  -----

InterventionManager::InterventionManager (const scnXml::Interventions& intervElt, Population& population) :
    nextTimed(0)
{
    if( intervElt.getVaccine().present() ){
        const scnXml::Vaccine& vacc = intervElt.getVaccine().get();
        if( vacc.getContinuous().size() + vacc.getTimed().size() > 0 ){
            activeInterventions.set (Interventions::VACCINE, true);
            // read descriptions:
            Host::Vaccine::init( vacc );
            // continuous deployments:
            //FIXME
            // timed deployments:
            typedef scnXml::Vaccine::TimedSequence::const_iterator It;
            for( It it = vacc.getTimed().begin(); it != vacc.getTimed().end(); ++it ){
                timed.push_back( new TimedMassCumIntervention( *it, population, &Host::Human::massVaccinate, &Host::Human::hasVaccineProtection ) );
            }
        }
    }
    if( intervElt.getMDA().present() ){
        const scnXml::MDA& mda = intervElt.getMDA().get();
        if( mda.getTimed().size() > 0 ){
            activeInterventions.set( Interventions::MDA, true );
            // read description:
            if( TimeStep::interval == 5 ){
                if( mda.getDescription().present() ){
                    cerr << "warning: MDA description not expected for 5-day timestep" << endl;
                }
            }else{
                if( !mda.getDescription().present() ){
                    throw util::xml_scenario_error( "error: MDA description required for 1-day timestep" );
                }
                Clinical::ESCaseManagement::initMDA( mda.getDescription().get() );
            }
            // timed deployments:
            typedef scnXml::MDA::TimedSequence::const_iterator It;
            for( It it = mda.getTimed().begin(); it != mda.getTimed().end(); ++it ){
                timed.push_back( new TimedMassIntervention( *it, population, &Host::Human::massDrugAdministration ) );
            }
        }
    }
    if( intervElt.getIPT().present() ){
        const scnXml::IPT& ipt = intervElt.getIPT().get();
        // read description (note: this is required by IPT WIH model, so do it even if there is no deployment)
        WithinHost::DescriptiveIPTWithinHost::init( ipt.getDescription() );
        // continuous deployments:
        //FIXME
        // timed deployments:
        typedef scnXml::IPT::TimedSequence::const_iterator It;
        for( It it = ipt.getTimed().begin(); it != ipt.getTimed().end(); ++it ){
            timed.push_back( new TimedMassCumIntervention( *it, population, &Host::Human::IPTiTreatment, &Host::Human::hasIPTiProtection ) );
        }
    }
    /*FIXME
    if (scenario->getInterventions().getContinuous().present()) {
        const scnXml::ContinuousInterv& contI = scenario->getInterventions().getContinuous().get();
        if (contI.getITN().size())
            activeInterventions.set (Interventions::ITN, true);
        if (contI.getIpti().size())
            activeInterventions.set (Interventions::IPTI, true);
        if (contI.getCohort().size())
            activeInterventions.set (Interventions::COHORT, true);
    }

    if (scenario->getInterventions().getTimed().present()) {
        const scnXml::Timed::InterventionSequence& interventionSeq =
            scenario->getInterventions().getTimed().get().getIntervention();
        for (scnXml::Timed::InterventionConstIterator it (interventionSeq.begin()); it != interventionSeq.end(); ++it) {
            TimeStep time = TimeStep(it->getTime());
            if (timedInterventions.count (time)) {
                ostringstream msg;
                msg << "Error: multiple timed interventions with time: " << time;
                throw util::xml_scenario_error (msg.str());
            }
            timedInterventions[time] = & (*it);

            if (it->getChangeHS().present())
                activeInterventions.set (Interventions::CHANGE_HS, true);
            if (it->getChangeEIR().present())
                activeInterventions.set (Interventions::CHANGE_EIR, true);
            if (it->getIpti().present())
                activeInterventions.set (Interventions::IPTI, true);
            if (it->getITN().present())
                activeInterventions.set (Interventions::ITN, true);
            if (it->getIRS().present())
                activeInterventions.set (Interventions::IRS, true);
            if (it->getVectorAvailability().present())
                activeInterventions.set (Interventions::VEC_AVAIL, true);
            if (it->getImmuneSuppression().present())
                activeInterventions.set (Interventions::IMMUNE_SUPPRESSION, true);
            if (it->getCohort().present())
                activeInterventions.set (Interventions::COHORT, true);
            if (it->getLarviciding().present())
                activeInterventions.set (Interventions::LARVICIDING, true);
            if (it->getInsertR_0Case().present()){
                activeInterventions.set (Interventions::R_0_CASE, true);
                // uses vaccines but see note in Vaccine::initParameters()
                // activeInterventions.set (Interventions::VACCINE, true);
            }
            if (it->getImportedInfectionsPerThousandHosts().present())
                activeInterventions.set (Interventions::IMPORTED_INFECTIONS, true);
            if (it->getUninfectVectors().present())
                activeInterventions.set (Interventions::UNINFECT_VECTORS, true);

        }
    }
    */
    // make sure the list ends with something always in the future, so we don't
    // have to check nextTimed is within range:
    timed.push_back( new DummyIntervention() );
    
    //FIXME: sort lists, earliest to latest
}

void InterventionManager::deploy() {
    // deploy continuous interventions (FIXME)
    // deploy timed interventions
    while( timed[nextTimed].time <= TimeStep::simulation ){
        timed[nextTimed].deploy();
        nextTimed += 1;
    }
}

}
