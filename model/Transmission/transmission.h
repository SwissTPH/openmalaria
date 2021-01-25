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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
 * USA.
 */

#ifndef Hmod_Transmission
#define Hmod_Transmission

#include "Transmission/NonVectorModel.h"
#include "Transmission/VectorModel.h"
#include "Transmission/Anopheles/AnophelesModel.h"
#include "Transmission/Anopheles/SimpleMPDAnophelesModel.h"

namespace OM
{
namespace Transmission
{

// static void loadVectorInterventions(const scnXml::Interventions& intervElt, VectorModel &model)
// {
//     if (intervElt.getVectorPop().present())
//     {
//         typedef scnXml::VectorPop::InterventionSequence SeqT;
//         const SeqT &seq = intervElt.getVectorPop().get().getIntervention();
//         size_t instance = 0;
//         for (auto it = seq.begin(), end = seq.end(); it != end; ++it)
//         {
//             const scnXml::VectorIntervention &elt = *it;
//             if (elt.getTimed().present())
//             {
//                 transmission.initVectorInterv(elt.getDescription().getAnopheles(), instance, elt.getName());

//                 const scnXml::TimedBaseList::DeploySequence &seq = elt.getTimed().get().getDeploy();
//                 for (auto it = seq.begin(); it != seq.end(); ++it)
//                 {
//                     SimDate date = UnitParse::readDate(it->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
//                     timed.push_back(unique_ptr<TimedDeployment>(new TimedVectorDeployment(date, instance)));
//                 }
//                 instance++;
//             }
//         }
//     }
//     if (intervElt.getAddNonHumanHosts().present())
//     {
//         typedef scnXml::AddNonHumanHosts::NonHumanHostsSequence SeqT;
//         const SeqT &seq = intervElt.getAddNonHumanHosts().get().getNonHumanHosts();
//         size_t instance = 0;
//         for (auto it = seq.begin(), end = seq.end(); it != end; ++it)
//         {
//             const scnXml::NonHumanHosts2 &elt = *it;
//             if (elt.getTimed().present())
//             {
//                 transmission.initAddNonHumanHostsInterv(elt.getDescription().getAnopheles(), elt.getName());
//                 for (const scnXml::Deploy2 deploy : elt.getTimed().get().getDeploy())
//                 {
//                     SimDate date = UnitParse::readDate(deploy.getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
//                     SimTime lifespan = UnitParse::readDuration(deploy.getLifespan(), UnitParse::NONE);
//                     timed.push_back(unique_ptr<TimedDeployment>(new TimedAddNonHumanHostsDeployment(date, elt.getName(), lifespan)));
//                 }
//                 instance++;
//             }
//         }
//     }
//     if (intervElt.getNonHumanHostsModifications().present())
//     {
//         typedef scnXml::NonHumanHostsModifications::InterventionSequence SeqT;
//         const SeqT &seq = intervElt.getNonHumanHostsModifications().get().getIntervention();
//         size_t instance = 0;
//         for (auto it = seq.begin(), end = seq.end(); it != end; ++it)
//         {
//             const scnXml::NonHumanHostsIntervention &elt = *it;
//             if (elt.getTimed().present())
//             {
//                 const scnXml::DecayFunction &decay = elt.getDecay();
//                 transmission.initNonHumanHostsInterv(elt.getDescription().getAnopheles(), decay, instance, elt.getNonHumanHostsName());
//                 const scnXml::TimedBaseList::DeploySequence &seq = elt.getTimed().get().getDeploy();
//                 for (auto it = seq.begin(); it != seq.end(); ++it)
//                 {
//                     SimDate date = UnitParse::readDate(it->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
//                     timed.push_back(
//                         unique_ptr<TimedDeployment>(new TimedNonHumanHostsDeployment(date, instance, elt.getNonHumanHostsName())));
//                 }
//                 instance++;
//             }
//         }
//     }
//     if (intervElt.getVectorTrap().present())
//     {
//         size_t instance = 0;
//         for (const scnXml::VectorTrap &trap : intervElt.getVectorTrap().get().getIntervention())
//         {
//             transmission.initVectorTrap(trap.getDescription(), instance, trap.getName());
//             if (trap.getTimed().present())
//             {
//                 for (const scnXml::Deploy1 deploy : trap.getTimed().get().getDeploy())
//                 {
//                     SimDate date = UnitParse::readDate(deploy.getTime(), UnitParse::STEPS);
//                     double ratio = deploy.getRatioToHumans();
//                     SimTime lifespan = UnitParse::readDuration(deploy.getLifespan(), UnitParse::NONE);
//                     timed.push_back(unique_ptr<TimedDeployment>(new TimedTrapDeployment(date, instance, ratio, lifespan)));
//                 }
//             }
//             instance += 1;
//         }
//     }
// }

static SimulationMode readMode(const string &str)
{
    if (str == "forced")
        return forcedEIR;
    else if (str == "dynamic")
        return dynamicEIR;
    else
        // Note: originally 3 (transientEIRknown) could be specified; now it's
        // set automatically.
        throw util::xml_scenario_error(string("mode attribute invalid: ").append(str));
}

static NonVectorModel *createNonVectorModel(const scnXml::Entomology &entoData)
{
    const scnXml::NonVector &nonVectorData = entoData.getNonVector().get();

    int interventionMode = readMode(entoData.getMode());
    int eipDuration = nonVectorData.getEipDuration();

    vector<double> initialisationEIR;
    initialisationEIR.assign(sim::stepsPerYear(), 0.0);

    return new NonVectorModel(initialisationEIR, interventionMode, entoData, nonVectorData, eipDuration);
}

static bool anophelesCompare(const scnXml::AnophelesParams &a1, const scnXml::AnophelesParams &a2)
{
    return (a1.getSeasonality().getAnnualEIR().get() > a2.getSeasonality().getAnnualEIR().get());
}

static VectorModel *createVectorModel(const scnXml::Entomology &entoData, const scnXml::Interventions &intervElt, int populationSize)
{
    const scnXml::Vector &vectorData = entoData.getVector().get();
    scnXml::Vector::AnophelesSequence anophelesList = vectorData.getAnopheles();

    int interventionMode = readMode(entoData.getMode());

    vector<double> initialisationEIR;
    initialisationEIR.assign(sim::stepsPerYear(), 0.0);

    vector<Anopheles::AnophelesModel *> species;
    map<string, size_t> speciesIndex;

    size_t numSpecies = anophelesList.size();
    cout << "numSpecies " << numSpecies << endl;
    if (numSpecies < 1) throw util::xml_scenario_error("Can't use Vector model without data for at least one anopheles species!");

    // Sort Anopheles by Decreasing EIR
    sort(anophelesList.begin(), anophelesList.end(), anophelesCompare);

    PerHostAnophParams::initReserve(numSpecies);
    species.resize(numSpecies);

    for (size_t i = 0; i < numSpecies; ++i)
    {
        auto elt = anophelesList[i];

        Anopheles::AnophelesModel *anophModel;

        if (util::ModelOptions::option(util::VECTOR_LIFE_CYCLE_MODEL))
        {
            throw util::xml_scenario_error("VECTOR_LIFE_CYCLE_MODEL not yet "
                                           "implemented. Use VECTOR_SIMPLE_MPD_MODEL instead.");
            // TODO
            //  * Note: this model is older than SimpleMPD and more complicated.
            //  * Difficulties are in parameterisation and estimation of resources.
            // if (!lcOpt.present())
            //     throw util::xml_scenario_error(
            //         "VECTOR_LIFE_CYCLE_MODEL: requires <lifeCycle> element with "
            //         "model parameters for each anopheles species");
            // emergence = unique_ptr<EmergenceModel>( new LCEmergence() );
            // emergence->initLifeCycle( lcOpt.get() );
        }
        else if (util::ModelOptions::option(util::VECTOR_SIMPLE_MPD_MODEL))
        {
            if (!elt.getSimpleMPD().present())
                throw util::xml_scenario_error("VECTOR_SIMPLE_MPD_MODEL: requires <simpleMPD> element with "
                                               "model parameters for each anopheles species");

            const scnXml::SimpleMPD& smpd = elt.getSimpleMPD().get();

            SimTime developmentDuration = SimTime::fromDays(smpd.getDevelopmentDuration().getValue());
            if (!(developmentDuration > SimTime::zero()))
                throw util::xml_scenario_error("entomology.vector.simpleMPD.developmentDuration: "
                    "must be positive");
            double probPreadultSurvival = smpd.getDevelopmentSurvival().getValue();
            if (!(0.0 <= probPreadultSurvival && probPreadultSurvival <= 1.0))
                throw util::xml_scenario_error("entomology.vector.simpleMPD.developmentSurvival: "
                    "must be a probability (in range [0,1]");
            double fEggsLaidByOviposit = smpd.getFemaleEggsLaidByOviposit().getValue();
            if (!(fEggsLaidByOviposit > 0.0))
                throw util::xml_scenario_error("entomology.vector.simpleMPD.femaleEggsLaidByOviposit: "
                    "must be positive");

            anophModel = new Anopheles::SimpleMPDAnophelesModel(developmentDuration, probPreadultSurvival, fEggsLaidByOviposit);
        }
        else
            anophModel = new Anopheles::AnophelesModel();

        PerHostAnophParams::init(elt.getMosq());

        species[i] = anophModel;
        species[i]->initialise(i, elt, initialisationEIR, populationSize);

        string name = elt.getMosquito();
        speciesIndex[name] = i;
    }

    if (interventionMode == forcedEIR)
    {
        // We don't need these anymore (now we have initialisationEIR); free memory
        numSpecies = 0;
        species.clear();
        speciesIndex.clear();
    }

    return new VectorModel(initialisationEIR, interventionMode, species, speciesIndex, populationSize);
}

///@brief Creation, destruction and checkpointing
//@{
/// Creates a derived class
static TransmissionModel *createTransmissionModel(const scnXml::Entomology &entoData, const scnXml::Interventions &intervElt,
                                                  int populationSize)
{
    // Entomology contains either a list of at least one anopheles or a list of at
    // least one EIRDaily.
    TransmissionModel *model;
    if (entoData.getVector().present())
        model = createVectorModel(entoData, intervElt, populationSize);
    else
    {
        if (!entoData.getNonVector().present()) // should be a validation error, but anyway...
            throw util::xml_scenario_error("Neither vector nor non-vector data present in the XML!");
        if (util::ModelOptions::option(util::VECTOR_LIFE_CYCLE_MODEL) || util::ModelOptions::option(util::VECTOR_SIMPLE_MPD_MODEL))
            throw util::xml_scenario_error("VECTOR_*_MODEL is only compatible with the vector model (and non-vector data is present).");
        model = createNonVectorModel(entoData);
    }

    if (entoData.getScaledAnnualEIR().present())
    {
        model->scaleEIR(entoData.getScaledAnnualEIR().get() / model->annualEIR);
        assert(util::vectors::approxEqual(model->annualEIR, entoData.getScaledAnnualEIR().get()));
    }

    if (util::CommandLine::option(util::CommandLine::PRINT_ANNUAL_EIR))
    {
        // Note: after internal scaling (which doesn't imply exit)
        // but before external scaling.
        cout << "Total annual EIR: " << model->annualEIR << endl;
    }

    return model;
}
} // namespace Transmission
} // namespace OM

#endif