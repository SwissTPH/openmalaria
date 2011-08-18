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

#include "Transmission/VectorModel.h"
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

double VectorModel::invMeanPopAvail (const std::list<Host::Human>& population, int populationSize) {
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

void VectorModel::ctsCbN_v0 (ostream& stream) {
    for (size_t i = 0; i < numSpecies; ++i)
        stream << '\t' << species[i].getLastN_v0()/TimeStep::interval;
}
void VectorModel::ctsCbP_A (ostream& stream) {
    for (size_t i = 0; i < numSpecies; ++i)
        stream << '\t' << species[i].getLastVecStat(Vector::PA)/TimeStep::interval;
}
void VectorModel::ctsCbP_df (ostream& stream) {
    for (size_t i = 0; i < numSpecies; ++i)
        stream << '\t' << species[i].getLastVecStat(Vector::PDF)/TimeStep::interval;
}
void VectorModel::ctsCbP_dif (ostream& stream) {
    for (size_t i = 0; i < numSpecies; ++i)
        stream << '\t' << species[i].getLastVecStat(Vector::PDIF)/TimeStep::interval;
}
void VectorModel::ctsCbN_v (ostream& stream) {
    for (size_t i = 0; i < numSpecies; ++i)
        stream << '\t' << species[i].getLastVecStat(Vector::NV)/TimeStep::interval;
}
void VectorModel::ctsCbO_v (ostream& stream) {
    for (size_t i = 0; i < numSpecies; ++i)
        stream << '\t' << species[i].getLastVecStat(Vector::OV)/TimeStep::interval;
}
void VectorModel::ctsCbS_v (ostream& stream) {
    for (size_t i = 0; i < numSpecies; ++i)
        stream << '\t' << species[i].getLastVecStat(Vector::SV)/TimeStep::interval;
}
void VectorModel::ctsCbResAvailability (ostream& stream) {
    for (size_t i = 0; i < numSpecies; ++i)
        stream << '\t' << species[i].getResAvailability()/TimeStep::interval;
}
void VectorModel::ctsCbResRequirements (ostream& stream) {
    for (size_t i = 0; i < numSpecies; ++i)
        stream << '\t' << species[i].getResRequirements()/TimeStep::interval;
}

const string& reverseLookup (const map<string,size_t>& m, size_t i) {
    for ( map<string,size_t>::const_iterator it = m.begin(); it != m.end(); ++it ) {
        if ( it->second == i )
            return it->first;
    }
    throw TRACED_EXCEPTION_DEFAULT( "reverseLookup: key not found" );        // shouldn't ever happen
}
VectorModel::VectorModel (const scnXml::Vector vectorData, int populationSize)
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
    species.resize (numSpecies, Vector::SpeciesModel(&_ITNParams));

    for (size_t i = 0; i < numSpecies; ++i) {
        string name = species[i].initialise (anophelesList[i],
                                             initialisationEIR,
                                             nonHumanHostPopulations,
                                             populationSize);
        speciesIndex[name] = i;
    }

    // Calculate total annual EIR
    annualEIR = vectors::sum( initialisationEIR );


    if( interventionMode == forcedEIR ) {
        // We don't need these anymore (now we have initialisationEIR); free memory
        numSpecies = 0;
        species.clear();
        speciesIndex.clear();
    }


    // -----  Continuous reporting  -----
    ostringstream ctsNv0, ctsPA, ctsPdf, ctsPdif, ctsNv, ctsOv, ctsSv, ctsRA, ctsRR;
    // Output in order of species so that (1) we can just iterate through this
    // list when outputting and (2) output is in order specified in XML.
    for (size_t i = 0; i < numSpecies; ++i) {
        // Unfortunately, we need to reverse-lookup the name.
        const string& name = reverseLookup( speciesIndex, i );
        ctsNv0<<"\tN_v0("<<name<<")";
        ctsPA<<"\tP_A("<<name<<")";
        ctsPdf<<"\tP_df("<<name<<")";
        ctsPdif<<"\tP_dif("<<name<<")";
        ctsNv<<"\tN_v("<<name<<")";
        ctsOv<<"\tO_v("<<name<<")";
        ctsSv<<"\tS_v("<<name<<")";
        ctsRA<<"\tres avail("<<name<<")";
        ctsRR<<"\tres req("<<name<<")";
    }
    using Monitoring::Continuous;
    Continuous::registerCallback( "N_v0", ctsNv0.str(), MakeDelegate( this, &VectorModel::ctsCbN_v0 ) );
    Continuous::registerCallback( "P_A", ctsPA.str(), MakeDelegate( this, &VectorModel::ctsCbP_A ) );
    Continuous::registerCallback( "P_df", ctsPdf.str(), MakeDelegate( this, &VectorModel::ctsCbP_df ) );
    Continuous::registerCallback( "P_dif", ctsPdif.str(), MakeDelegate( this, &VectorModel::ctsCbP_dif ) );
    Continuous::registerCallback( "N_v", ctsNv.str(), MakeDelegate( this, &VectorModel::ctsCbN_v ) );
    Continuous::registerCallback( "O_v", ctsOv.str(), MakeDelegate( this, &VectorModel::ctsCbO_v ) );
    Continuous::registerCallback( "S_v", ctsSv.str(), MakeDelegate( this, &VectorModel::ctsCbS_v ) );
    Continuous::registerCallback( "resource availability", ctsRA.str(), MakeDelegate( this, &VectorModel::ctsCbResAvailability ) );
    Continuous::registerCallback( "resource requirements", ctsRR.str(), MakeDelegate( this, &VectorModel::ctsCbResRequirements ) );
}
VectorModel::~VectorModel () {
}

void VectorModel::setupNv0 (const std::list<Host::Human>& population, int populationSize) {
    double iMPA = invMeanPopAvail(population, populationSize);
    for (size_t i = 0; i < numSpecies; ++i) {
        species[i].setupNv0 (i, population, populationSize, iMPA);
    }
    simulationMode = forcedEIR;   // now we should be ready to start
}


void VectorModel::scaleEIR (double factor) {
    //FIXME: this needs revision; rename XML elt to "overrideAnnualEIR"? Store annual EIR per species.
    for ( size_t i = 0; i < numSpecies; ++i )
        species[i].scaleEIR( factor );
    vectors::scale( initialisationEIR, factor );
    annualEIR = vectors::sum( initialisationEIR );
}
void VectorModel::scaleXML_EIR (scnXml::EntoData& ed, double factor) const {
    //FIXME: this needs revision; in fact it's completely redundant now but maybe re can re-use code to save resource availability?
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


TimeStep VectorModel::minPreinitDuration () {
    if ( interventionMode == forcedEIR ) {
        return TimeStep(0);
    }
    // Data is summed over 5 years; add an extra 50 for stabilization.
    // 50 years seems a reasonable figure from a few tests
    return TimeStep::fromYears( 55 );
}
TimeStep VectorModel::expectedInitDuration (){
    return TimeStep::fromYears( 1 );
}

TimeStep VectorModel::initIterate () {
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
        throw TRACED_EXCEPTION("Transmission warmup exceeded 10 iterations!",util::Error::VectorWarmup);
    }
    
    // Time to let parameters settle after each iteration. I would expect one year
    // to be enough (but I may be wrong).
    if( needIterate ){
        return TimeStep::fromYears( 1 ) + TimeStep::fromYears(5);  // stabilization + 5 years data-collection time
    }
    return TimeStep::fromYears( 1 );
}

double VectorModel::calculateEIR(PerHost& host, double ageYears) {
    host.update(_ITNParams);
    if (simulationMode == forcedEIR){
        return initialisationEIR[TimeStep::simulation % TimeStep::stepsPerYear]
               * host.relativeAvailabilityHetAge (ageYears);
    }else{
        assert( simulationMode == dynamicEIR );
        double simEIR = 0.0;
        for (size_t i = 0; i < numSpecies; ++i) {
            simEIR += species[i].calculateEIR (i, host);
        }
        simEIR *= host.relativeAvailabilityAge (ageYears);
        return simEIR;
    }
}


// Every Global::interval days:
void VectorModel::vectorUpdate (const std::list<Host::Human>& population, int populationSize) {
    if( simulationMode == dynamicEIR ){
        double iMPA = invMeanPopAvail(population, populationSize);
        for (size_t i = 0; i < numSpecies; ++i){
            species[i].advancePeriod (population, populationSize, i, iMPA);
        }
    }
}
void VectorModel::update (const std::list<Host::Human>& population, int populationSize) {
    TransmissionModel::updateKappa( population );
}

inline void assertIsDynamic( int interventionMode ){
    if( interventionMode != dynamicEIR ){
        throw xml_scenario_error("vector interventions can only be used in dynamic transmission mode");
    }
}
void VectorModel::setITNDescription (const scnXml::ITNDescription& elt){
    assertIsDynamic( interventionMode );
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
void VectorModel::setIRSDescription (const scnXml::IRS& elt){
    assertIsDynamic( interventionMode );
    PerHost::setIRSDescription (elt);
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
void VectorModel::setVADescription (const scnXml::VectorDeterrent& elt){
    assertIsDynamic( interventionMode );
    PerHost::setVADescription (elt);
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

void VectorModel::intervLarviciding (const scnXml::Larviciding& anoph) {
    assertIsDynamic( interventionMode );
    assert(false);
    /*FIXME
    const scnXml::Larviciding::AnophelesSequence& seq = anoph.getAnopheles();
    for (scnXml::Larviciding::AnophelesSequence::const_iterator it = seq.begin(); it != seq.end(); ++it)
        species[getSpeciesIndex(it->getMosquito())].intervLarviciding(*it);
    */
}
void VectorModel::uninfectVectors() {
    for (size_t i = 0; i < numSpecies; ++i)
        species[i].uninfectVectors();
}

void VectorModel::summarize (Monitoring::Survey& survey) {
    TransmissionModel::summarize (survey);

    for (map<string,size_t>::const_iterator it = speciesIndex.begin(); it != speciesIndex.end(); ++it)
        species[it->second].summarize (it->first, survey);
}


void VectorModel::checkpoint (istream& stream) {
    TransmissionModel::checkpoint (stream);
    initIterations & stream;
    util::checkpoint::checkpoint (species, stream, Vector::SpeciesModel (&_ITNParams));
}
void VectorModel::checkpoint (ostream& stream) {
    TransmissionModel::checkpoint (stream);
    initIterations & stream;
    species & stream;
}

}
}
