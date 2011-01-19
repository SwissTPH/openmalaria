/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
 OpenMalaria is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
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

#include "WithinHost/DescriptiveIPTWithinHost.h"
#include "WithinHost/Infection/DescriptiveIPTInfection.h"
#include "util/random.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "inputData.h"
#include "Monitoring/Surveys.h"
#include "PopulationStats.h"

#include <cmath>

namespace OM { namespace WithinHost {
    using namespace util;

bool DescriptiveIPTWithinHost::iptActive = false;

// -----  static data  -----

DescriptiveIPTWithinHost::IPTiEffects DescriptiveIPTWithinHost::iptiEffect;


// -----  init  -----

void DescriptiveIPTWithinHost::init () {
  iptActive = util::ModelOptions::option( IPTI_SP_MODEL );
  if (!iptActive) {
    if (InputData.getActiveInterventions()[Interventions::IPTI])
      throw util::xml_scenario_error ("IPTI interventions require IPT_SP_MODEL option");
    return;
  }
  const scnXml::Interventions& xmlInterventions = InputData().getInterventions();
  if (!xmlInterventions.getDescriptions().getIptiDescription().present()) {
    throw util::xml_scenario_error ("IPT_SP_MODEL requires iptiDescription");
  }
  
  if (TimeStep::interval != 5)
    throw util::xml_scenario_error ("IPT code only supports using an interval of 5");
  if (util::ModelOptions::option (util::INCLUDES_PK_PD)) {
      throw util::xml_scenario_error ("DescriptiveIPTWithinHost not intended to work with DrugAction");
      // The IPT code has its own implementation of non-instantaneous drug action (SPAction, etc).
  }
  
  // --- IptiDescription begin ---
  const scnXml::IptDescription& xmlIPTI = xmlInterventions.getDescriptions().getIptiDescription().get();
  
  iptiEffect = static_cast<IPTiEffects>(xmlIPTI.getIptiEffect());
  // --- IptiDescription end ---
  
  DescriptiveIPTInfection::initParameters(xmlInterventions);
}

void DescriptiveIPTWithinHost::cleanup () {
  if (!iptActive) return;
  DescriptiveIPTInfection::cleanup();
}

DescriptiveIPTWithinHost::DescriptiveIPTWithinHost () :
    _SPattenuationt(TimeStep::never), _lastSPDose (TimeStep::never), _lastIptiOrPlacebo (TimeStep::never),
    _cumulativeInfections(0)
{
}


// -----  Simple infection adders/removers  -----

void DescriptiveIPTWithinHost::newInfection(){
    ++PopulationStats::totalInfections;
  if (_MOI < MAX_INFECTIONS) {
    _cumulativeInfections++;
    infections.push_back(new DescriptiveIPTInfection(_lastSPDose));
    _MOI++;
    ++PopulationStats::allowedInfections;
  }
  assert( _MOI == infections.size() );
}
void DescriptiveIPTWithinHost::loadInfection(istream& stream){
    infections.push_back(new DescriptiveIPTInfection(stream));
}

// -----  Events setting _lastSPDose  -----

void DescriptiveIPTWithinHost::clearInfections (bool isSevere) {
  TimeStep fortnight = TimeStep::fromDaysNearest(14.0);	// round to nearest
  if (isSevere) {
  } else if(TimeStep::simulation-_lastIptiOrPlacebo <= fortnight) {
          // IPTi trials used quinine for fevers within 14 days of an ipti or placebo dose   
  } else if(TimeStep::simulation-_lastSPDose <=  fortnight) {
          /*
    second line used if fever within 14 days of SP dose (ipti or treatment)
    
    if this code is to survive, then the iptiEffect values should be 
    symbolic constants
    however: code is dead (only used for repeat experiments) anyway
          */
  } else if( iptiEffect ==  PLACEBO_SP ||  iptiEffect ==  IPT_SP) {
      // add 1 to delay start until next update since effect of SP already happened this time step
    _lastSPDose=TimeStep::simulation+TimeStep(1);
  } else if( iptiEffect ==  PLACEBO_CLEAR_INFECTIONS ||  iptiEffect ==  IPT_CLEAR_INFECTIONS) {
  } else if(iptiEffect >= IPT_SEASONAL_MIN && iptiEffect < IPT_MAX) {
  } else {
      // add 1 to delay start until next update since effect of SP already happened this time step
    _lastSPDose=TimeStep::simulation+TimeStep(1);
    // SPAction will first act at the beginning of the next Global::interval
  }
  clearAllInfections();
}

void DescriptiveIPTWithinHost::deployIptDose (Monitoring::AgeGroup ageGroup, bool inCohort) {
  // assumes 5-day intervals and Niakhar seasonality
  // These numbers, should have MAX = MIN + 18 (modulo 73).
  static int IPT_MIN_INTERVAL[9] = { 43, 49, 55, 61, 67, 37, 31, 25, 19 };
  static int IPT_MAX_INTERVAL[9] = { 61, 67,  0,  6, 12, 55, 49, 43, 31 };
  
  if (iptiEffect >= IPT_SEASONAL_MIN && iptiEffect <= IPT_SEASONAL_MAX) {
    int yearInterval = TimeStep::simulation % TimeStep::stepsPerYear;
    int min = IPT_MIN_INTERVAL[iptiEffect-14];
    int max = IPT_MAX_INTERVAL[iptiEffect-14];
    // We're using modular arithmatic here, to represent a time period 5*18 days long.
    if (min < max) {
      if (yearInterval < min || yearInterval >= max)
	return;
    } else {
      if (yearInterval < min && yearInterval >= max)
	return;
    }
  }
  
    _lastIptiOrPlacebo=TimeStep::simulation;
    /*
    iptiEffect denotes treatment or placebo group
    and also the treatment given when sick (trial-dependent)
    */
    if (iptiEffect >= IPT_MIN) {
	_lastSPDose=TimeStep::simulation;
	Monitoring::Surveys.getSurvey(inCohort).reportIPTDoses (ageGroup, 1);
    }
}

void DescriptiveIPTWithinHost::IPTiTreatment (Monitoring::AgeGroup ageGroup, bool inCohort) {
  //Set the last SP Dose given for the eligible humans - is this all we need to do?
  
  _lastIptiOrPlacebo = TimeStep::simulation;
  
  // iptiEffect denotes treatment or placebo group
  // and also the treatment given when sick (trial-dependent)
  if (iptiEffect >= IPT_MIN){
    _lastSPDose = TimeStep::simulation;
    Monitoring::Surveys.getSurvey(inCohort).reportIPTDoses (ageGroup, 1);
  }
}
bool DescriptiveIPTWithinHost::hasIPTiProtection (TimeStep maxInterventionAge) const{
    return _lastIptiOrPlacebo + maxInterventionAge > TimeStep::simulation;
}


// -----  update (from SPAction)  -----

bool DescriptiveIPTWithinHost::eventSPClears (DescriptiveInfection* inf){
  DescriptiveIPTInfection* iptInf = dynamic_cast<DescriptiveIPTInfection*> (inf);
  assert (iptInf);	// code error if this is anything else
  return iptInf->eventSPClears(_lastSPDose);
}

// -----  density calculation  -----

void DescriptiveIPTWithinHost::IPTattenuateAsexualMinTotalDensity () {
  //Note: the _cumulativeInfections>0 check is probably unintended, but was extracted from other logic and put here to preserve results.
  if (util::ModelOptions::option (util::ATTENUATION_ASEXUAL_DENSITY) && _cumulativeInfections > 0) {
    if (_SPattenuationt > TimeStep::simulation && totalDensity < 10) {
      totalDensity = 10;
      _cumulativeY += 10;
    }
  }
}

void DescriptiveIPTWithinHost::IPTattenuateAsexualDensity (DescriptiveInfection* inf) {
  if (!(util::ModelOptions::option (util::ATTENUATION_ASEXUAL_DENSITY))) return;
  
  DescriptiveIPTInfection* iptInf = dynamic_cast<DescriptiveIPTInfection*> (inf);
  assert (iptInf);	// code error if this is anything else
  if (iptInf->doSPAttenuation()) {
    double attFact = iptInf->asexualAttenuation();
    timeStepMaxDensity *= attFact;
    _SPattenuationt = std::max(_SPattenuationt, iptInf->getAsexualAttenuationEndDate());
  }
}


// -----  Data checkpointing  -----

void DescriptiveIPTWithinHost::checkpoint (istream& stream) {
    DescriptiveWithinHostModel::checkpoint (stream);
    _SPattenuationt & stream;
    _lastSPDose & stream; 
    _lastIptiOrPlacebo & stream;
    _cumulativeInfections & stream;
}
void DescriptiveIPTWithinHost::checkpoint (ostream& stream) {
    DescriptiveWithinHostModel::checkpoint (stream);
    _SPattenuationt & stream;
    _lastSPDose & stream; 
    _lastIptiOrPlacebo & stream;
    _cumulativeInfections & stream;
}

} }