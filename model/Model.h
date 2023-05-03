#ifndef MODEL_H
#define MODEL_H

#include "Global.h"

#include "Population.h"
#include "Transmission/transmission.h"

namespace OM {

struct Model
{
    Model(Population *population, Transmission::TransmissionModel *transmission, bool surveyOnlyNewEp, SimTime humanWarmupLength) : 
	    population(unique_ptr<Population>(population)),
	    transmission(unique_ptr<Transmission::TransmissionModel>(transmission)),
	    surveyOnlyNewEp(surveyOnlyNewEp),
	    humanWarmupLength(humanWarmupLength)
	    {}

    std::unique_ptr<Population> population;
    std::unique_ptr<Transmission::TransmissionModel> transmission;

    bool surveyOnlyNewEp;
    SimTime humanWarmupLength;
};

namespace model
{
	Model* create(const scnXml::Scenario &scenario);

	// Internal simulation loop
	void run(Model &model, SimTime &endTime, SimTime &estEndTime, string phase);
}

}

#endif 