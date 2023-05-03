#include "Global.h"

#include "Model.h"

#include "Parameters.h"

#include "mon/Continuous.h"
#include "mon/management.h"

#include "Host/NeonatalMortality.h"
#include "Host/Human.h"
#include "Host/InfectionIncidenceModel.h"
#include "Host/WithinHost/WHInterface.h"

#include "interventions/InterventionManager.hpp"

#include "Clinical/ClinicalModel.h"

namespace OM {

using mon::Continuous;
using interventions::InterventionManager;
using Transmission::TransmissionModel;

namespace model
{
Model* create(const scnXml::Scenario &scenario)
{
    sim::init(scenario);
    
    // 1) elements with no dependencies on other elements initialised here:
    // sim::init( *scenario );  // also reads survey dates
    Parameters parameters( scenario.getModel().getParameters() );     // depends on nothing
    WithinHost::Genotypes::init( scenario );
    
    util::master_RNG.seed( scenario.getModel().getParameters().getIseed(), 0 ); // Init RNG with Iseed
    util::ModelOptions::init( scenario.getModel().getModelOptions() );
    
    // 2) elements depending on only elements initialised in (1):
    WithinHost::diagnostics::init( parameters, scenario ); // Depends on Parameters
    mon::initReporting( scenario ); // Reporting init depends on diagnostics and monitoring
    
    // Init models used by humans
    Transmission::PerHost::init( scenario.getModel().getHuman().getAvailabilityToMosquitoes() );
    Host::InfectionIncidenceModel::init( parameters );
    WithinHost::WHInterface::init( parameters, scenario );
    Clinical::ClinicalModel::init( parameters, scenario );
    Host::NeonatalMortality::init( scenario.getModel().getClinical() );
    AgeStructure::init( scenario.getDemography() );

    // 3) elements depending on other elements; dependencies on (1) are not mentioned:
    // Transmission model initialisation depends on Transmission::PerHost and
    // genotypes (both from Human, from Population::init()) and
    // mon::AgeGroup (from Surveys.init()):
    // Note: PerHost dependency can be postponed; it is only used to set adultAge
    size_t popSize = scenario.getDemography().getPopSize();
    Population *population = population::createPopulation(popSize);
    TransmissionModel *transmission = Transmission::createTransmissionModel(scenario.getEntomology(), popSize);

    // Depends on transmission model (for species indexes):
    // MDA1D may depend on health system (too complex to verify)
    interventions::InterventionManager::init(scenario.getInterventions(), *population, *transmission );
    Clinical::ClinicalModel::setHS( scenario.getHealthSystem() ); // Depends on interventions, PK/PD (from humanPop)
    mon::initCohorts( scenario.getMonitoring() ); // Depends on interventions

    bool surveyOnlyNewEp = scenario.getMonitoring().getSurveyOptions().getOnlyNewEpisode();

    sim::s_t0 = sim::zero();
    sim::s_t1 = sim::zero();
    
    // Make sure warmup period is at least as long as a human lifespan, as the
    // length required by vector warmup, and is a whole number of years.
    SimTime humanWarmupLength = sim::maxHumanAge();
    if( humanWarmupLength < transmission->minPreinitDuration() ){
        cerr << "Warning: human life-span (" << sim::inYears(humanWarmupLength);
        cerr << ") shorter than length of warm-up requested by" << endl;
        cerr << "transmission model (" << sim::inYears(transmission->minPreinitDuration());
        cerr << "). Transmission may be unstable; perhaps use forced" << endl;
        cerr << "transmission (mode=\"forced\") or a longer life-span." << endl;
        humanWarmupLength = transmission->minPreinitDuration();
    }
    humanWarmupLength = sim::fromYearsI( static_cast<int>(ceil(sim::inYears(humanWarmupLength))) );

    return new Model(population, transmission, surveyOnlyNewEp, humanWarmupLength);
}

void print_progress(int lastPercent, SimTime &estEndTime)
{
    int percent = (sim::now() * 100) / estEndTime;
    if( percent != lastPercent ) {   // avoid huge amounts of output for performance/log-file size reasons
        lastPercent = percent;
        cerr << "\r" << percent << "%\t" << flush;
    }
}

void print_errno()
{
    if( errno != 0 )
    {
       char err[256];
       sprintf(err, "t = %d Please report! Error: ", int(sim::now()));
       std::perror(err);
       errno = 0;
    }
}

// Internal simulation loop
void run(Model &model, SimTime &endTime, SimTime &estEndTime, string phase)
{
    static int lastPercent = -1;

    Population &population = *model.population;
    TransmissionModel &transmission = *model.transmission;

    if (util::CommandLine::option(util::CommandLine::VERBOSE)) cout << "Starting " << phase << "..." << endl;

    while (sim::now() < endTime)
    {
        if (util::CommandLine::option(util::CommandLine::VERBOSE) && sim::intervDate() > 0)
            cout << "Time step: " << sim::now() / sim::oneTS() << ", internal days: " << sim::now() << " | " << estEndTime << ", Intervention Date: " << sim::intervDate() << endl;

        // Monitoring. sim::now() gives time of end of last step,
        // and is when reporting happens in our time-series.
        Continuous.update( population );
        if( sim::intervDate() == mon::nextSurveyDate() ){
            for(Host::Human &human : population.humans)
                Host::human::summarize(human, model.surveyOnlyNewEp);
            transmission.summarize();
            mon::concludeSurvey();
        }
        
        // Deploy interventions, at time sim::now().
        InterventionManager::deploy( population.humans, transmission );
        
        // Time step updates. Time steps are mid-day to mid-day.
        // sim::ts0() gives the date at the start of the step, sim::ts1() the date at the end.
        sim::start_update();

        // This should be called before humans contract new infections in the simulation step.
        // This needs the whole population (it is an approximation before all humans are updated).
        transmission.vectorUpdate(population.humans);
        
        // NOTE: no neonatal mortalities will occur in the first 20 years of warmup
        // (until humans old enough to be pregnate get updated and can be infected).
        Host::NeonatalMortality::update (population.humans);
        
        for (Host::Human& human : population.humans)
        {
            if (human.dateOfBirth + sim::maxHumanAge() >= model.humanWarmupLength) // this is last time of possible update
                Host::human::update(human, transmission);
        }
       
        population.regularize();
        
        // Doesn't matter whether non-updated humans are included (value isn't used
        // before all humans are updated).
        transmission.updateKappa(population.humans);

        sim::end_update();

        print_progress(lastPercent, estEndTime);
        print_errno();
    }

    if (util::CommandLine::option(util::CommandLine::VERBOSE)) cout << "Finishing " << phase << "..." << endl;
}

}

}
