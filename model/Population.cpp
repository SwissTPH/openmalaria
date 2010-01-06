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
#include "util/gsl.h"
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

#include <cmath>

namespace OM
{

// -----  Population: static data / methods  -----

int Population::workUnitIdentifier;

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
    
    workUnitIdentifier = InputData.get_wu_id();
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
    
    workUnitIdentifier & stream;
    
    // Check scenario.xml and checkpoint files correspond:
    if (workUnitIdentifier != InputData.get_wu_id())
	throw util::checkpoint_error ("invalid work-unit identifier");
}
void Population::staticCheckpoint (ostream& stream)
{
    Host::NeonatalMortality::staticCheckpoint (stream);
    Clinical::ClinicalModel::staticCheckpoint (stream);
    PkPd::PkPdModel::staticCheckpoint (stream);
    
    workUnitIdentifier & stream;
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
    //TODO: cleanup & check (this code duplicates old behaviour)
    for (int iage = AgeStructure::getMaxTimestepsPerLife() - 2; iage >= 0; iage--) {
	int targetPop = AgeStructure::targetCumPop (2+iage, populationSize);
	while (cumulativePop < targetPop) {
	    newHuman (-iage);
	    ++cumulativePop;
	}
    }
    /*    for (int j = 1;j < maxTimestepsPerLife; j++) {
    int iage = maxTimestepsPerLife - j - 1;
    int targetPop = (int) floor (cumAgeProp[j] * populationSize + 0.5);
    while (cumulativePop < targetPop) {
   if (InitPopOpt && iage > 0) {} // only those with age 0 should be created here
       else newHuman (-iage);
       ++cumulativePop;
       }
       }*/
    
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
        //FIXME: The +2 here is to replicate old results. I think it's wrong though. Also, it looks
	// like this code assumes the maximum age of indivs is maxTimestepsPerLife not
	// Global::maxAgeIntervals.
        if (cumPop > AgeStructure::targetCumPop (age + 2, targetPop)) {
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
        massIntervention (interv->getMDA().get(), &Host::Human::clearInfections);
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


// -----  Population::AgeStructure  -----

int Population::AgeStructure::maxTimestepsPerLife;

double Population::AgeStructure::ageGroupBounds[ngroups+1];
double Population::AgeStructure::ageGroupPercent[ngroups];

double Population::AgeStructure::M1[ngroups];
double Population::AgeStructure::M2[ngroups];
double Population::AgeStructure::M[ngroups];
double Population::AgeStructure::pred[ngroups];

double Population::AgeStructure::mu0;
double Population::AgeStructure::mu1;
double Population::AgeStructure::alpha0;
double Population::AgeStructure::alpha1;
double Population::AgeStructure::rho;

vector<double> Population::AgeStructure::cumAgeProp;


void Population::AgeStructure::init () {
    maxTimestepsPerLife = maxLifetimeDays / Global::interval;
    cumAgeProp.resize (maxTimestepsPerLife);
    
    estimateRemovalRates();
    calcCumAgeProp();
}

int Population::AgeStructure::targetCumPop (int ageTSteps, int targetPop)
{
    return (int) floor (cumAgeProp[maxTimestepsPerLife+1-ageTSteps] * targetPop + 0.5);
}

void Population::AgeStructure::estimateRemovalRates ()
{
    // mu1, alpha1: These are estimated here
    // mu0, alpha0: These are fixed alpha0=4.0, mu0 calculated from the other parameters
    // rho - population growth rate (input)
    
    //Get lower and upper age bounds for age groups and cumulative precentage of population from field data
    double sumperc = 0.0;
    const scnXml::AgeGroupPerC::GroupSequence& group = InputData.getDemography().getAgeGroup().getGroup();
    if (group.size() < ngroups - 1) {
	ostringstream msg;
	msg << "expected " << ngroups - 1 << " elements of \"group\" in demography->ageGroup (in scenario.xml)";
	throw util::xml_scenario_error (msg.str());
    }
    //Add age group for first month of life
    ageGroupBounds[0] = 0.0;
    ageGroupBounds[1] = 1.0 / 12.0;
    ageGroupPercent[0] = 0.0;
    for (int i = 1;i < ngroups; i++) {
	ageGroupBounds[i+1] = group[i-1].getUpperbound();
	ageGroupPercent[i] = group[i-1].getPoppercent();
	sumperc += ageGroupPercent[i];
    }
    sumperc = 100.0 / sumperc; // multiplier to get percentages
    for (int i = 0;i < ngroups; i++) {
	ageGroupPercent[i]  = ageGroupPercent[i] * sumperc;
    }
    /*
    RSS between observed and predicted log percentage of population in age groups
    is minimised for values of mu1 and alpha1
    calls setDemoParameters to calculate the RSS
    */
    /* NOTE: unused --- why?
    double tol = 0.00000000001;
    int maxcal = 500000;
    int npar = 2;
    int iw = 3; */
    double p1 = 0.371626412;
    double p2 = 0.841209593;
    // returns "double rss":
    gsl::minimizeCalc_rss (&setDemoParameters, p1, p2);
}

// Static method used by estimateRemovalRates
double Population::AgeStructure::setDemoParameters (double param1, double param2)
{
    rho = InputData.get_growthrate() * (0.01 * Global::yearsPerInterval);
    if (rho != 0.0)
	// Issue: in this case the total population size differs from populationSize,
	// however, some code currently uses this as the total population size.
	throw util::xml_scenario_error ("Population growth rate provided.");
    
    const double IMR = 0.1;
    double M_inf = -log (1 - IMR);
    
    mu1 = exp (param1) / 100;
    alpha1 = exp (param2) / 100;
    alpha0 = 4.0;
    mu0 = (M_inf - mu1 * (exp (alpha1 * 0.5) - 1) * alpha0) /
    (alpha1 * (1 - exp (-alpha0 * 0.5)));
    
    double sumpred = 0.0;
    for (int i = 0; i < ngroups - 1; i++) {
	double midpt = (ageGroupBounds[i+1] + ageGroupBounds[i]) * 0.5;
	M1[i] = mu0 * (1.0 - exp (-alpha0 * midpt)) / alpha0;
	M2[i] = mu1 * (exp (alpha1 * midpt) - 1.0) / alpha1;
	M[i]  = M1[i] + M2[i];
	pred[i] = (ageGroupBounds[i+1] - ageGroupBounds[i])
	* exp (-rho * midpt - M[i]);
	sumpred += pred[i];
    }
    for (int i = 0; i < ngroups - 1; i++) {
	pred[i] = pred[i] / sumpred * 100.0;
    }
    double L_inf = exp (-rho * 0.5 - M[1]);
    double M_nn  = -log (1.0 - 0.4 * (1 - exp (-M[1])));
    double L1    = 1.0 / 12.0 * exp (-rho / 24.0 - M_nn);
    double perc_inf = ageGroupPercent[0] + ageGroupPercent[1];
    ageGroupPercent[0] = perc_inf * L1 / L_inf;
    ageGroupPercent[1] = perc_inf - ageGroupPercent[0];
    
    double valsetDemoParameters = 0.0;
    for (int i = 0; i < ngroups - 1; i++) {
	double residual = log (pred[i]) - log (ageGroupPercent[i]);
	valsetDemoParameters += residual * residual;
    }
    return valsetDemoParameters;
}

void Population::AgeStructure::calcCumAgeProp ()
{
    cumAgeProp[0] = 0.0;
    for (int j = 1;j < maxTimestepsPerLife; j++) {
	double ageYears = (maxTimestepsPerLife - j - 1) * Global::yearsPerInterval;
	double M1s = (mu0 * (1.0 - exp (-alpha0 * ageYears)) / alpha0);
	double M2s = (mu1 * (exp (alpha1 * ageYears) - 1.0) / alpha1);
	double Ms = M1s + M2s;
	double predperc = exp (-rho * ageYears - Ms);
	if (j < maxTimestepsPerLife - Global::maxAgeIntervals) {
	    predperc = 0.0;
	}
	cumAgeProp[j] = cumAgeProp[j-1] + predperc;
    }
    double totalCumPC = cumAgeProp[maxTimestepsPerLife-1];
    for (int j = 1;j < maxTimestepsPerLife; j++) {
	//Scale using the total cumAgeProp
	cumAgeProp[j] = cumAgeProp[j] / totalCumPC;
    }
}

}
