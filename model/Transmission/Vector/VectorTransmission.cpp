/* This file is part of OpenMalaria.
 *
 * Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#include "Transmission/Vector/VectorTransmission.h"
#include "inputData.h"
#include "Monitoring/Continuous.h"
#include "util/vectors.h"
#include "util/ModelOptions.h"

#include <fstream>
#include <map>
#include <cmath>

namespace OM {
namespace Transmission {
using namespace OM::util;

double VectorTransmission::invMeanPopAvail (const std::list<Host::Human>& population, int populationSize) {
    double sumRelativeAvailability = 0.0;
    for (std::list<Host::Human>::const_iterator h = population.begin(); h != population.end(); ++h){
        sumRelativeAvailability += h->perHostTransmission.relativeAvailabilityAge (h->getAgeInYears());
    }
    if( sumRelativeAvailability > 0.0 ){
        return populationSize / sumRelativeAvailability;     // 1 / mean-rel-avail
    }else{
        // value should be unimportant when no humans are available, though inf/nan is not acceptable
        return 1.0;
    }
}

void VectorTransmission::ctsCbN_v0 (ostream& stream) {
    for (size_t i = 0; i < numSpecies; ++i)
        stream << '\t' << species[i].getLastN_v0()/TimeStep::interval;
}
void VectorTransmission::ctsCbN_v (ostream& stream) {
    for (size_t i = 0; i < numSpecies; ++i)
        stream << '\t' << species[i].getLastN_v()/TimeStep::interval;
}
void VectorTransmission::ctsCbO_v (ostream& stream) {
    for (size_t i = 0; i < numSpecies; ++i)
        stream << '\t' << species[i].getLastO_v()/TimeStep::interval;
}
void VectorTransmission::ctsCbS_v (ostream& stream) {
    for (size_t i = 0; i < numSpecies; ++i)
        stream << '\t' << species[i].getLastS_v()/TimeStep::interval;
}
const string& reverseLookup (const map<string,size_t>& m, size_t i) {
    for ( map<string,size_t>::const_iterator it = m.begin(); it != m.end(); ++it ) {
        if ( it->second == i )
            return it->first;
    }
    throw util::traced_exception( "reverseLookup: key not found" );        // shouldn't ever happen
}

VectorTransmission::VectorTransmission (const scnXml::Vector vectorData, int populationSize)
    : initIterations(0), numSpecies(0)
{
    // Each item in the AnophelesSequence represents an anopheles species.
    // TransmissionModel::createTransmissionModel checks length of list >= 1.
    const scnXml::Vector::AnophelesSequence anophelesList = vectorData.getAnopheles();
    const scnXml::Vector::NonHumanHostsSequence nonHumansList = vectorData.getNonHumanHosts();

    map<string, double> nonHumanHostPopulations;

    for (size_t i = 0; i<nonHumansList.size(); i++)
        nonHumanHostPopulations[nonHumansList[i].getName()] = nonHumansList[i].getNumber();

    numSpecies = anophelesList.size();
    if (numSpecies < 1)
        throw util::xml_scenario_error ("Can't use Vector model without data for at least one anopheles species!");
    species.resize (numSpecies, VectorAnopheles(&_ITNParams));

    for (size_t i = 0; i < numSpecies; ++i) {
        string name = species[i].initialise (anophelesList[i],
                                             initialisationEIR,
                                             nonHumanHostPopulations,
                                             populationSize);
        speciesIndex[name] = i;
    }

    // Calculate total annual EIR
    annualEIR = vectors::sum( initialisationEIR );


    if( interventionMode == equilibriumMode ) {
        // We don't need these anymore (now we have initialisationEIR); free memory
        numSpecies = 0;
        species.clear();
        speciesIndex.clear();
    }


    // -----  Continuous reporting  -----
    ostringstream ctsNv0, ctsNv, ctsOv, ctsSv;
    // Output in order of species so that (1) we can just iterate through this
    // list when outputting and (2) output is in order specified in XML.
    for (size_t i = 0; i < numSpecies; ++i) {
        // Unfortunately, we need to reverse-lookup the name.
        const string& name = reverseLookup( speciesIndex, i );
        ctsNv0<<"\tN_v0("<<name<<")";
        ctsNv<<"\tN_v("<<name<<")";
        ctsOv<<"\tO_v("<<name<<")";
        ctsSv<<"\tS_v("<<name<<")";
    }
    using Monitoring::Continuous;
    Continuous::registerCallback( "N_v0", ctsNv0.str(), MakeDelegate( this, &VectorTransmission::ctsCbN_v0 ) );
    Continuous::registerCallback( "N_v", ctsNv.str(), MakeDelegate( this, &VectorTransmission::ctsCbN_v ) );
    Continuous::registerCallback( "O_v", ctsOv.str(), MakeDelegate( this, &VectorTransmission::ctsCbO_v ) );
    Continuous::registerCallback( "S_v", ctsSv.str(), MakeDelegate( this, &VectorTransmission::ctsCbS_v ) );
}
VectorTransmission::~VectorTransmission () {
}

void VectorTransmission::setupNv0 (const std::list<Host::Human>& population, int populationSize) {
    double iMPA = invMeanPopAvail(population, populationSize);
    for (size_t i = 0; i < numSpecies; ++i) {
        species[i].setupNv0 (i, population, populationSize, iMPA);
    }
}


void VectorTransmission::scaleEIR (double factor) {
    for ( size_t i = 0; i < numSpecies; ++i )
        species[i].scaleEIR( factor );
    vectors::scale( initialisationEIR, factor );
    annualEIR = vectors::sum( initialisationEIR );
}
void VectorTransmission::scaleXML_EIR (scnXml::EntoData& ed, double factor) const {
    // XML values are exponentiated; so we add some factor to existing a0 values:
    double add_to_a0 = std::log( factor );

    assert( ed.getVector().present() );
    scnXml::Vector::AnophelesSequence& anophelesList = ed.getVector().get().getAnopheles();

    for ( scnXml::Vector::AnophelesIterator it = anophelesList.begin();
            it != anophelesList.end(); ++it ) {
        if( it->getEIR().present() ){
            double old_a0 = it->getEIR().get().getA0();
            it->getEIR().get().setA0( old_a0 + add_to_a0 );
        }else{
            // schema should enforce that one of the two is here
            assert( it->getMonthlyEIR().present() );
            double oldAnnual = it->getMonthlyEIR().get().getAnnualEIR();
            it->getMonthlyEIR().get().setAnnualEIR( oldAnnual * factor );
        }
    }
}


TimeStep VectorTransmission::minPreinitDuration () {
    if ( interventionMode == equilibriumMode ) {
        return TimeStep(0);
    }
    // Data is summed over 5 years; add an extra 50 for stabilization.
    // 50 years seems a reasonable figure from a few tests
    return TimeStep::fromYears( 55 );
}
TimeStep VectorTransmission::expectedInitDuration (){
    return TimeStep::fromYears( 1 );
}

TimeStep VectorTransmission::initIterate () {
    if( interventionMode != dynamicEIR ) {
        // allow forcing equilibrium mode like with non-vector model
        return TimeStep(0); // no initialization to do
    }
    if( initIterations < 0 ){
        assert( interventionMode = dynamicEIR );
        simulationMode = dynamicEIR;
        return TimeStep(0);
    }
    
    ++initIterations;
    
    bool needIterate = false;
    for (size_t i = 0; i < numSpecies; ++i) {
        needIterate = needIterate || species[i].vectorInitIterate ();
    }
    if( needIterate == false ){
        initIterations = -1;
    }
    if( initIterations > 10 ){
        throw util::traced_exception("Transmission warmup exceeded 10 iterations!",util::Error::VectorWarmup);
    }
    
    // Time to let parameters settle after each iteration. I would expect one year
    // to be enough (but I may be wrong).
    if( needIterate ){
        return TimeStep::fromYears( 1 ) + TimeStep::fromYears(5);  // stabilization + 5 years data-collection time
    }
    return TimeStep::fromYears( 1 );
}

double VectorTransmission::calculateEIR(PerHostTransmission& host, double ageYears) {
    host.update(_ITNParams);
    if (simulationMode == equilibriumMode){
        return initialisationEIR[TimeStep::simulation % TimeStep::stepsPerYear]
               * host.relativeAvailabilityHetAge (ageYears);
    }else{      // dynamicMode
        double simEIR = 0.0;
        for (size_t i = 0; i < numSpecies; ++i) {
            simEIR += species[i].calculateEIR (i, host);
        }
        simEIR *= host.relativeAvailabilityAge (ageYears);
        return simEIR;
    }
}


// Every Global::interval days:
void VectorTransmission::vectorUpdate (const std::list<Host::Human>& population, int populationSize) {
    double iMPA = invMeanPopAvail(population, populationSize);
    for (size_t i = 0; i < numSpecies; ++i){
        species[i].advancePeriod (population, populationSize, i, simulationMode == dynamicEIR, iMPA);
    }
}
void VectorTransmission::update (const std::list<Host::Human>& population, int populationSize) {
    TransmissionModel::updateKappa( population );
}

void VectorTransmission::setITNDescription (const scnXml::ITNDescription& elt){
    if( interventionMode == equilibriumMode ){
        throw xml_scenario_error("vector interventions can only be used in dynamic transmission mode (mode=4)");
    }
    double proportionUse = _ITNParams.init( elt );
    typedef scnXml::ITNDescription::AnophelesParamsSequence AP;
    const AP& ap = elt.getAnophelesParams();
    if( ap.size() != numSpecies ){
        throw util::xml_scenario_error(
            "ITNDescription.anophelesParams: must have one element for each mosquito species described in entomology"
        );
    }
    for( AP::const_iterator it = ap.begin(); it != ap.end(); ++it ){
        species[getSpeciesIndex(it->getMosquito())].setITNDescription (_ITNParams, *it, proportionUse);
    }
}
void VectorTransmission::setIRSDescription (const scnXml::IRS& elt){
    if( interventionMode == equilibriumMode ){
        throw xml_scenario_error("vector interventions can only be used in dynamic transmission mode (mode=4)");
    }
    PerHostTransmission::setIRSDescription (elt);
    typedef scnXml::IRS::AnophelesParamsSequence AP;
    const AP& ap = elt.getAnophelesParams();
    if( ap.size() != numSpecies ){
        throw util::xml_scenario_error(
            "IRS.anophelesParams: must have one element for each mosquito species described in entomology"
        );
    }
    for( AP::const_iterator it = ap.begin(); it != ap.end(); ++it ) {
        species[getSpeciesIndex(it->getMosquito())].setIRSDescription (*it);
    }
}
void VectorTransmission::setVADescription (const scnXml::VectorDeterrent& elt){
    if( interventionMode == equilibriumMode ){
        throw xml_scenario_error("vector interventions can only be used in dynamic transmission mode (mode=4)");
    }
    PerHostTransmission::setVADescription (elt);
    typedef scnXml::VectorDeterrent::AnophelesParamsSequence AP;
    const AP& ap = elt.getAnophelesParams();
    if( ap.size() != numSpecies ){
        throw util::xml_scenario_error(
            "vectorDeterrent.anophelesParams: must have one element for each mosquito species described in entomology"
        );
    }
    for( AP::const_iterator it = ap.begin(); it != ap.end(); ++it ) {
        species[getSpeciesIndex(it->getMosquito())].setVADescription (*it);
    }
}

void VectorTransmission::intervLarviciding (const scnXml::Larviciding& anoph) {
    if( interventionMode == equilibriumMode ){
        throw xml_scenario_error("vector interventions can only be used in dynamic transmission mode (mode=4)");
    }
    assert(false);
    /*FIXME
    const scnXml::Larviciding::AnophelesSequence& seq = anoph.getAnopheles();
    for (scnXml::Larviciding::AnophelesSequence::const_iterator it = seq.begin(); it != seq.end(); ++it)
        species[getSpeciesIndex(it->getMosquito())].intervLarviciding(*it);
    */
}
void VectorTransmission::uninfectVectors() {
    for (size_t i = 0; i < numSpecies; ++i)
        species[i].uninfectVectors();
}

void VectorTransmission::summarize (Monitoring::Survey& survey) {
    TransmissionModel::summarize (survey);

    for (map<string,size_t>::const_iterator it = speciesIndex.begin(); it != speciesIndex.end(); ++it)
        species[it->second].summarize (it->first, survey);
}


void VectorTransmission::checkpoint (istream& stream) {
    TransmissionModel::checkpoint (stream);
    initIterations & stream;
    util::checkpoint::checkpoint (species, stream, VectorAnopheles (&_ITNParams));
}
void VectorTransmission::checkpoint (ostream& stream) {
    TransmissionModel::checkpoint (stream);
    initIterations & stream;
    species & stream;
}

}
}
