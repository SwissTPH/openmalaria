/*

 This file is part of OpenMalaria.

 Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine

 OpenMalaria is free software; you can redistribute it and/or modify
 it under the terms of the GNu General Public License as published by
 the Free Software Foundation; either version 2 of the License, or (at
 your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

*/
#include "Simulation.h"
#include "Population.h"
#include "inputData.h"
#include "Surveys.h"

#include "Transmission/TransmissionModel.h"

#include "Host/Human.h"
#include "Host/NeonatalMortality.h"
#include "Clinical/ImmediateOutcomes.h"
#include "Pathogenesis/PathogenesisModel.h"
#include "PkPd/PkPdModel.h"

#include "util/errors.hpp"
#include "util/ModelOptions.hpp"

namespace OM
{

// -----  Population: static data / methods  -----

#ifdef OMP_CSV_REPORTING
ofstream csvReporting;
#endif


void Population::init()
{
    Host::Human::initHumanParameters();
    Host::NeonatalMortality::init();
    PkPd::PkPdModel::init();
#ifdef OMP_CSV_REPORTING
    csvReporting.open ("population.csv", ios::app);
#endif
    
    AgeStructure::init ();
}

void Population::clear()
{
    PkPd::PkPdModel::cleanup ();
    Host::Human::clear();
#ifdef OMP_CSV_REPORTING
    csvReporting.close();
#endif
}

void Population::staticCheckpoint (istream& stream)
{
    Host::NeonatalMortality::staticCheckpoint (stream);
    Clinical::ClinicalModel::staticCheckpoint (stream);
    PkPd::PkPdModel::staticCheckpoint (stream);
}
void Population::staticCheckpoint (ostream& stream)
{
    Host::NeonatalMortality::staticCheckpoint (stream);
    Clinical::ClinicalModel::staticCheckpoint (stream);
    PkPd::PkPdModel::staticCheckpoint (stream);
}


// -----  non-static methods: creation/destruction, checkpointing  -----

Population::Population()
        : populationSize (InputData.get_populationsize())
{
    _transmissionModel = Transmission::TransmissionModel::createTransmissionModel();
}

Population::~Population()
{
    for (HumanIter iter = population.begin(); iter != population.end(); ++iter) {
        iter->destroy();
    }
    delete _transmissionModel;
}

void Population::checkpoint (istream& stream)
{
    size_t popSize; // must be type of population.size()
    popSize & stream;
    if (popSize != size_t (populationSize))
        throw util::checkpoint_error ("population size exceeds that given in scenario.xml");
    while (popSize > 0 && !stream.eof()) {
        // Note: calling this constructor of Host::Human is slightly wasteful, but avoids the need for another
        // ctor and leaves less opportunity for uninitialized memory.
        population.push_back (Host::Human (*_transmissionModel, 0, 0));
        population.back() & stream;
        --popSize;
    }
    if (int (population.size()) != populationSize)
        throw util::checkpoint_error ("can't read whole population (out of data)");
}
void Population::checkpoint (ostream& stream)
{
    population.size() & stream;
    for (HumanIter iter = population.begin(); iter != population.end(); ++iter)
        (*iter) & stream;
}

void Population::preMainSimInit ()
{
    _transmissionModel->initMainSimulation();
    Clinical::ClinicalModel::initMainSimulation();
}

void Population::createInitialHumans ()
{
    int cumulativePop = 0;
    for (int iage = AgeStructure::getMaxTimestepsPerLife() - 1; iage >= 0; iage--) {
	int targetPop = AgeStructure::targetCumPop (iage, populationSize);
	while (cumulativePop < targetPop) {
	    newHuman (-iage);
	    ++cumulativePop;
	}
    }
    
    // Vector setup dependant on human population
    _transmissionModel->updateAgeCorrectionFactor (population, populationSize);
    _transmissionModel->setupNv0 (population, populationSize);
}


// -----  non-static methods: simulation loop  -----

void Population::newHuman (int dob)
{
    population.push_back (Host::Human (*_transmissionModel, dob, Global::simulationTime));
}

void Population::update1()
{
    _transmissionModel->updateAgeCorrectionFactor (population, populationSize);

    Host::NeonatalMortality::update (population);
    // This should be called before humans contract new infections in the simulation step.
    _transmissionModel->vectorUpdate (population, Global::simulationTime);

    //NOTE: other parts of code are not set up to handle changing population size. Also
    // populationSize is assumed to be the _actual and exact_ population size by other code.
    //targetPop is the population size at time t allowing population growth
    //int targetPop = (int) (populationSize * exp (AgeStructure::rho * Global::simulationTime));
    int targetPop = populationSize;
    int cumPop = 0;

    // Update each human in turn
    //std::cout<<" time " <<t<<std::endl;
    HumanIter last = population.end();
    --last;
    for (HumanIter iter = population.begin(); iter != population.end();) {
        // Update human, and remove if too old:
        if (iter->update (Global::simulationTime, _transmissionModel)) {
            iter->destroy();
            iter = population.erase (iter);
            continue;
        }

        //BEGIN Population size & age structure
        ++cumPop;
        int age = (Global::simulationTime - iter->getDateOfBirth());

        // if (Actual number of people so far > target population size for this age) ...
        if (cumPop > AgeStructure::targetCumPop (age, targetPop)) {
            --cumPop;
            iter->destroy();
            iter = population.erase (iter);
            continue;
        }
        //END Population size & age structure
        ++iter;
    } // end of per-human updates

    // increase population size to targetPop
    while (cumPop < targetPop) {
        newHuman (Global::simulationTime);
        //++nCounter;
        ++cumPop;
    }

    _transmissionModel->updateKappa (population, Global::simulationTime);

#ifdef OMP_CSV_REPORTING
    if (Global::simulationTime % (Global::intervalsPerYear*5) == 0) {
        csvReporting << Global::simulationTime << ',';
        list<Host::Human>::reverse_iterator it = population.rbegin();
        for (double ageLim = 0; ageLim <= maxLifetimeDays / 365.0; ageLim += 1) {
            int counter = 0;
            while (it != population.rend() && it->getAgeInYears() < ageLim) {
                ++counter;
                ++it;
            }
            csvReporting << counter << ',';
        }
        csvReporting << endl;
    }
#endif
}


// -----  non-static methods: summarising and interventions  -----

void Population::newSurvey ()
{
    Survey& current = *Surveys.current;
    for (HumanIter iter = population.begin(); iter != population.end(); iter++) {
        iter->summarize (current);
    }
    _transmissionModel->summarize (current);
}

void Population::implementIntervention (int time)
{
    const scnXml::Intervention* interv = InputData.getInterventionByTime (time);
    if (interv == NULL)
        return;

    // Given an intervention descriptor for this time point, check which
    // interventions are included. All should check when data loading that they
    // have data if used according to getActiveInterventions().

    if (interv->getChangeHS().present()) {
        if (util::ModelOptions::option (util::CLINICAL_EVENT_SCHEDULER))
            throw util::xml_scenario_error ("Only ClinicalImmediateOutcomes is compatible with change of health-system intervention.");
        Clinical::OldCaseManagement::setHealthSystem(time);
    }

    if (interv->getChangeEIR().present()) {
        _transmissionModel->changeEIRIntervention (interv->getChangeEIR().get());
    }

    if (interv->getVaccinate().present()) {
        massIntervention (interv->getVaccinate().get(), &Host::Human::massVaccinate);
    }
    if (interv->getMDA().present()) {
	/* TODO: here we assume a 100% clearance rate for the MDA drug we use.
	This is not consistent with the way we treat according to the Health
	system description. The default clearance rate for MDA should be 100%
	since this simulates what was meant to happen in Garki.  We can change
	this by introducing an optional clearance rate that can be < 100%. */
	//TODO: MDA with the old treatment system should do this; with the drug model it needs some drugs to medicate.
	massIntervention (interv->getMDA().get().getCoverage(), &Host::Human::massDrugAdministration);
    }
    if (interv->getIpti().present()) {
        massIntervention (interv->getIpti().get(), &Host::Human::IPTiTreatment);
    }

    if (interv->getITN().present()) {
        massIntervention (interv->getITN().get(), &Host::Human::setupITN);
    }
    if (interv->getIRS().present()) {
        massIntervention (interv->getIRS().get(), &Host::Human::setupIRS);
    }
    if (interv->getVectorAvailability().present()) {
        massIntervention (interv->getVectorAvailability().get(), &Host::Human::setupVA);
    }

    if (interv->getLarviciding().present()) {
        _transmissionModel->intervLarviciding (interv->getLarviciding().get());
    }
}

void Population::massIntervention (const scnXml::Mass& mass, void (Host::Human::*intervention) ())
{
    double minAge = mass.getMinAge().present() ?
                    mass.getMinAge().get() : 0.0;
    double maxAge = mass.getMaxAge().present() ?
                    mass.getMaxAge().get() : 100.0;
    double coverage = mass.getCoverage();

    for (HumanIter iter = population.begin(); iter != population.end(); ++iter) {
        double ageYears = iter->getAgeInYears();
        if ( (ageYears > minAge) && (ageYears < maxAge) && gsl::rngUniform() < coverage)
            // This is UGLY syntax. It just means call intervention() on the human pointed by iter.
            ( (*iter).*intervention) ();
    }
}

}
