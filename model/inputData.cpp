/* This file is part of OpenMalaria.
 *
 * Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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


#include "inputData.h"
#include "scenario.hxx"
#include "util/BoincWrapper.h"
#include "util/errors.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <map>

namespace OM {
    /// Current schema version.
const int SCHEMA_VERSION = 12;
/** Oldest which current code is potentially compatible with
 * (provided the scenario.xml file references this version and doesn't use
 * members changed in newer versions). */
const int OLDEST_COMPATIBLE = 5;

// Initialization functions:

void InputDataType::initParameterValues()
{
    // initialize all to zero
    for (size_t i = 0; i < Params::MAX; ++i)
        parameterValues[i] = 0;
    // set parameters
    const scnXml::Parameters::ParameterSequence& paramSeq = parameters->getParameter();
    for (scnXml::Parameters::ParameterConstIterator it = paramSeq.begin(); it != paramSeq.end(); ++it) {
        int i = it->getNumber();
        if (i < 0 || i >= Params::MAX)
            std::cerr << "Warning: parameter with invalid index; ignoring." << std::endl;
        else
            parameterValues[i] = it->getValue();
    }
}

void InputDataType::initTimedInterventions()
{
    if (interventions->getContinuous().present()) {
        const scnXml::Continuous& contI = interventions->getContinuous().get();
        if (contI.getVaccine().size())
            activeInterventions.set (Interventions::VACCINE, true);
	if (contI.getITN().size())
	    activeInterventions.set (Interventions::ITN, true);
	if (contI.getIpti().size())
	    activeInterventions.set (Interventions::IPTI, true);
    }

    if (interventions->getTimed().present()) {
        const scnXml::Timed::InterventionSequence& interventionSeq =
            interventions->getTimed().get().getIntervention();
        for (scnXml::Timed::InterventionConstIterator it (interventionSeq.begin()); it != interventionSeq.end(); ++it) {
            int time = it->getTime();
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
            if (it->getVaccinate().present())
                activeInterventions.set (Interventions::VACCINE, true);
            if (it->getMDA().present())
                activeInterventions.set (Interventions::MDA, true);
            if (it->getIpti().present())
                activeInterventions.set (Interventions::IPTI, true);
            if (it->getITN().present())
                activeInterventions.set (Interventions::ITN, true);
            if (it->getIRS().present())
                activeInterventions.set (Interventions::IRS, true);
            if (it->getVectorAvailability().present())
                activeInterventions.set (Interventions::VEC_AVAIL, true);
            if (it->getLarviciding().present())
                activeInterventions.set (Interventions::LARVICIDING, true);
        }
    }
}


void InputDataType::createDocument (std::string lXmlFile)
{
    xmlFileName = lXmlFile;
    //Parses the document
    //NOTE: it'd be nice if this used Global::lookupResource for the schema.
# ifdef WITHOUT_BOINC
    // We don't need a checksum when not run with BOINC, so open by filename,
    // which allows using the filename to help find the schema file.
    // Note that the schema location can be set manually by passing properties,
    // but we won't necessarily have the right schema version associated with
    // the XML file in that case.
    scenario = (scnXml::parseScenario (lXmlFile.c_str())).release();
# else
    ifstream fileStream (lXmlFile.c_str());
    scenario = (scnXml::parseScenario (fileStream)).release();
    BoincWrapper::generateChecksum (fileStream);
    fileStream.close ();
#endif
    if (scenario->getSchemaVersion() < OLDEST_COMPATIBLE) {
        ostringstream msg;
        msg << "Input scenario.xml uses an outdated schema version; please update with SchemaTranslator. Current version: " << SCHEMA_VERSION;
        throw util::xml_scenario_error (msg.str());
    }
    if (scenario->getSchemaVersion() > SCHEMA_VERSION)
        throw util::xml_scenario_error ("Error: new schema version unsupported");

    monitoring = &scenario->getMonitoring();
    interventions = &scenario->getInterventions();
    entoData = &scenario->getEntoData();
    demography = &scenario->getDemography();
    if (scenario->getHealthSystem().present())
        healthSystem = &scenario->getHealthSystem().get();
    else
        healthSystem = NULL;
    eventScheduler = scenario->getEventScheduler().present() ?
                      &scenario->getEventScheduler().get() : NULL;
    parameters = &scenario->getParameters();
    //proteome = &scenario->getProteome();

    initParameterValues();
    initTimedInterventions();
}

void InputDataType::cleanDocument()
{
    if (documentChanged) {
        // get the "basename" (file name without path) of xmlFileName as a C string:
        const char* lastFS = strrchr (xmlFileName.c_str(), '/');
        const char* lastBS = strrchr (xmlFileName.c_str(), '\\');
        const char* baseName = lastBS > lastFS ? lastBS : lastFS;
        if (baseName == NULL) // no path separator found; use whole string
            baseName = xmlFileName.c_str();
        else
            ++baseName;  // start at next character

        ofstream outStream (baseName);
        ostringstream schema;
        schema << "scenario_" << SCHEMA_VERSION << ".xsd";

        xml_schema::NamespaceInfomap map;
        map[""].name = "";
        map[""].schema = schema.str();
        serializeScenario (outStream, *scenario, map);

        outStream.close();
    }

    // Destructors should handle cleanup
    if (scenario != NULL)
        delete scenario;
}


const scnXml::Scenario& InputDataType::getScenario() {
    return *scenario;
}

const scnXml::Monitoring& InputDataType::getMonitoring()
{
    return *monitoring;
}
const scnXml::Interventions& InputDataType::getInterventions()
{
    return *interventions;
}
const scnXml::EntoData& InputDataType::getEntoData()
{
    return *entoData;
}
const scnXml::Demography& InputDataType::getDemography()
{
    return *demography;
}
const scnXml::EventScheduler& InputDataType::getEventScheduler ()
{
    if (eventScheduler == NULL)
	throw util::xml_scenario_error ("EventScheduler requested but not in XML");
    return *eventScheduler;
}
const scnXml::HealthSystem& InputDataType::getHealthSystem()
{
    if (healthSystem == NULL)
        throw util::xml_scenario_error ("heathSystem element requested but not present");
    return *healthSystem;
}

scnXml::Scenario& InputDataType::getMutableScenario()
{
    return *scenario;
}

void InputDataType::changeHealthSystem (const scnXml::HealthSystem* hs)
{
    healthSystem = hs;
}

double InputDataType::getParameter (size_t i)
{
    return parameterValues[i];
}


// ----- Member access functions (bridges) -----
// This is largely unmodified from the old xerces version.

double InputDataType::get_detectionlimit()
{
    return  monitoring->getSurveys().getDetectionLimit();
}

int InputDataType::get_summary_option()
{
    return monitoring->getSurveys().getSummaryOption();
}

int InputDataType::get_mode()
{
    return scenario->getMode();
}

int InputDataType::get_assim_mode()
{
    return scenario->getAssimMode();
}

int InputDataType::get_wu_id()
{
    return scenario->getWuID();
}

double InputDataType::get_maximum_ageyrs()
{
    return scenario->getMaximumAgeYrs();
}

double InputDataType::get_lowerbound()
{
    return monitoring->getAgeGroup().getLowerbound();
}

const scnXml::Intervention* InputDataType::getInterventionByTime (int time)
{
    std::map<int, const scnXml::Intervention*>::iterator i = timedInterventions.find (time);
    if (i != timedInterventions.end())
        return i->second;
    else
        return NULL;
}
const bitset<Interventions::SIZE> InputDataType::getActiveInterventions ()
{
    return activeInterventions;
}

int InputDataType::get_analysis_no()
{
    return scenario->getAnalysisNo();
}

int InputDataType::get_populationsize()
{
    return scenario->getPopSize();
}


double InputDataType::get_growthrate()
{
    if (demography->getGrowthRate().present())
        return demography->getGrowthRate().get();
    else
        return 0.;
}

int InputDataType::get_latentp()
{
    return parameters->getLatentp();
}

int InputDataType::get_interval()
{
    return parameters->getInterval();
}

int InputDataType::getISeed()
{
    return parameters->getIseed();
}

InputDataType InputData;
}