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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "Transmission/VectorModel.h"
#include "Population.h"
#include "Host/Human.h"
#include "WithinHost/WHInterface.h"
#include "WithinHost/Genotypes.h"
#include "Monitoring/Continuous.h"
#include "util/vectors.h"
#include "util/ModelOptions.h"
#include "util/SpeciesIndexChecker.h"

#include <fstream>
#include <map>
#include <cmath>
#include <set>

namespace OM {
namespace Transmission {
using namespace OM::util;

void VectorModel::ctsCbN_v0 (ostream& stream) {
    for(size_t i = 0; i < numSpecies; ++i)
        stream << '\t' << species[i].getLastN_v0();
}
void VectorModel::ctsCbP_A (ostream& stream) {
    for(size_t i = 0; i < numSpecies; ++i)
        stream << '\t' << species[i].getLastVecStat(Anopheles::PA);
}
void VectorModel::ctsCbP_df (ostream& stream) {
    for(size_t i = 0; i < numSpecies; ++i)
        stream << '\t' << species[i].getLastVecStat(Anopheles::PDF);
}
void VectorModel::ctsCbP_dif (ostream& stream) {
    for(size_t i = 0; i < numSpecies; ++i)
        stream << '\t' << species[i].getLastVecStat(Anopheles::PDIF);
}
void VectorModel::ctsCbN_v (ostream& stream) {
    for(size_t i = 0; i < numSpecies; ++i)
        stream << '\t' << species[i].getLastVecStat(Anopheles::NV);
}
void VectorModel::ctsCbO_v (ostream& stream) {
    for(size_t i = 0; i < numSpecies; ++i)
        stream << '\t' << species[i].getLastVecStat(Anopheles::OV);
}
void VectorModel::ctsCbS_v (ostream& stream) {
    for(size_t i = 0; i < numSpecies; ++i)
        stream << '\t' << species[i].getLastVecStat(Anopheles::SV);
}
void VectorModel::ctsCbAlpha (const Population& population, ostream& stream){
    for( size_t i = 0; i < numSpecies; ++i){
        const Anopheles::PerHostBase& params = species[i].getHumanBaseParams();
        double total = 0.0;
        for(Population::ConstIter iter = population.cbegin(); iter != population.cend(); ++iter) {
            total += iter->perHostTransmission.entoAvailabilityFull( params, i, iter->age(sim::now()).inYears() );
        }
        stream << '\t' << total / population.size();
    }
}
void VectorModel::ctsCbP_B (const Population& population, ostream& stream){
    for( size_t i = 0; i < numSpecies; ++i){
	const Anopheles::PerHostBase& params = species[i].getHumanBaseParams();
        double total = 0.0;
        for(Population::ConstIter iter = population.cbegin(); iter != population.cend(); ++iter) {
            total += iter->perHostTransmission.probMosqBiting( params, i );
        }
        stream << '\t' << total / population.size();
    }
}
void VectorModel::ctsCbP_CD (const Population& population, ostream& stream){
    for( size_t i = 0; i < numSpecies; ++i){
	const Anopheles::PerHostBase& params = species[i].getHumanBaseParams();
        double total = 0.0;
        for(Population::ConstIter iter = population.cbegin(); iter != population.cend(); ++iter) {
            total += iter->perHostTransmission.probMosqResting( params, i );
        }
        stream << '\t' << total / population.size();
    }
}
void VectorModel::ctsNetInsecticideContent (const Population& population, ostream& stream){
//     double meanVar = 0.0;
//     int n = 0;
//     for(Population::ConstIter iter = population.cbegin(); iter != population.cend(); ++iter) {
//         if( iter->perHostTransmission.getITN().timeOfDeployment() >= SimTime::zero() ){
//             ++n;
//             meanVar += iter->perHostTransmission.getITN().getInsecticideContent(_ITNParams);
//         }
//     }
//     stream << '\t' << meanVar/n;
}
void VectorModel::ctsIRSInsecticideContent (const Population& population, ostream& stream) {
    //TODO(monitoring): work out how this applies when multiple IRS effects are allowed
//     double totalInsecticide = 0.0;
//     for(Population::ConstIter iter = population.cbegin(); iter != population.cend(); ++iter) {
//         totalInsecticide += iter->perHostTransmission.getIRS().getInsecticideContent(_IRSParams);
//     }
//     stream << '\t' << totalInsecticide / population.size();
}
void VectorModel::ctsIRSEffects (const Population& population, ostream& stream) {
    //TODO(monitoring): work out how this applies when multiple IRS effects are allowed
//     for( size_t i = 0; i < numSpecies; ++i ){
//         const interventions::IRSAnophelesParams& params = species[i].getHumanBaseParams().irs;
//         double totalRA = 0.0, totalPrePSF = 0.0, totalPostPSF = 0.0;
//         for(Population::ConstIter iter = population.cbegin(); iter != population.cend(); ++iter) {
//             totalRA += iter->perHostTransmission.getIRS().relativeAttractiveness(params);
//             totalPrePSF += iter->perHostTransmission.getIRS().preprandialSurvivalFactor(params);
//             totalPostPSF += iter->perHostTransmission.getIRS().postprandialSurvivalFactor(params);
//         }
//         stream << '\t' << totalRA / population.size()
//             << '\t' << totalPrePSF / population.size()
//             << '\t' << totalPostPSF / population.size();
//     }
}

void VectorModel::ctsCbResAvailability (ostream& stream) {
    for(size_t i = 0; i < numSpecies; ++i)
        stream << '\t' << species[i].getResAvailability();
}
void VectorModel::ctsCbResRequirements (ostream& stream) {
    for(size_t i = 0; i < numSpecies; ++i)
        stream << '\t' << species[i].getResRequirements();
}

const string& reverseLookup (const map<string,size_t>& m, size_t i) {
    for( map<string,size_t>::const_iterator it = m.begin(); it != m.end(); ++it ) {
        if ( it->second == i )
            return it->first;
    }
    throw TRACED_EXCEPTION_DEFAULT( "reverseLookup: key not found" );        // shouldn't ever happen
}

VectorModel::VectorModel (const scnXml::Entomology& entoData,
                          const scnXml::Vector vectorData, int populationSize) :
    TransmissionModel( entoData, WithinHost::Genotypes::N() ),
    initIterations(0), numSpecies(0)
{
    // Each item in the AnophelesSequence represents an anopheles species.
    // TransmissionModel::createTransmissionModel checks length of list >= 1.
    const scnXml::Vector::AnophelesSequence anophelesList = vectorData.getAnopheles();
    const scnXml::Vector::NonHumanHostsSequence nonHumansList = vectorData.getNonHumanHosts();

    map<string, double> nonHumanHostPopulations;

    for(size_t i = 0; i<nonHumansList.size(); i++)
        nonHumanHostPopulations[nonHumansList[i].getName()] = nonHumansList[i].getNumber();

    numSpecies = anophelesList.size();
    if (numSpecies < 1)
        throw util::xml_scenario_error ("Can't use Vector model without data for at least one anopheles species!");
    species.resize (numSpecies, AnophelesModel());

    for(size_t i = 0; i < numSpecies; ++i) {
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
    ostringstream ctsNv0, ctsPA, ctsPdf, ctsPdif,
        ctsNv, ctsOv, ctsSv,
        ctsAlpha, ctsPB, ctsPCD,
        ctsIRSEffects,
        ctsRA, ctsRR;
    // Output in order of species so that (1) we can just iterate through this
    // list when outputting and (2) output is in order specified in XML.
    for(size_t i = 0; i < numSpecies; ++i) {
        // Unfortunately, we need to reverse-lookup the name.
        const string& name = reverseLookup( speciesIndex, i );
        ctsNv0<<"\tN_v0("<<name<<")";
        ctsPA<<"\tP_A("<<name<<")";
        ctsPdf<<"\tP_df("<<name<<")";
        ctsPdif<<"\tP_dif("<<name<<")";
        ctsNv<<"\tN_v("<<name<<")";
        ctsOv<<"\tO_v("<<name<<")";
        ctsSv<<"\tS_v("<<name<<")";
        ctsAlpha<<"\talpha_i("<<name<<")";
        ctsPB<<"\tP_B("<<name<<")";
        ctsPCD<<"\tP_C*P_D("<<name<<")";
        ctsIRSEffects<<"\tIRS rel attr ("<<name<<")"
            <<"\tIRS preprand surv factor ("<<name<<")"
            <<"\tIRS postprand surv factor ("<<name<<")";
        ctsRA<<"\tres avail("<<name<<")";
        ctsRR<<"\tres req("<<name<<")";
    }
    using Monitoring::Continuous;
    Continuous.registerCallback( "N_v0", ctsNv0.str(), MakeDelegate( this, &VectorModel::ctsCbN_v0 ) );
    Continuous.registerCallback( "P_A", ctsPA.str(), MakeDelegate( this, &VectorModel::ctsCbP_A ) );
    Continuous.registerCallback( "P_df", ctsPdf.str(), MakeDelegate( this, &VectorModel::ctsCbP_df ) );
    Continuous.registerCallback( "P_dif", ctsPdif.str(), MakeDelegate( this, &VectorModel::ctsCbP_dif ) );
    Continuous.registerCallback( "N_v", ctsNv.str(), MakeDelegate( this, &VectorModel::ctsCbN_v ) );
    Continuous.registerCallback( "O_v", ctsOv.str(), MakeDelegate( this, &VectorModel::ctsCbO_v ) );
    Continuous.registerCallback( "S_v", ctsSv.str(), MakeDelegate( this, &VectorModel::ctsCbS_v ) );
    // availability to mosquitoes relative to other humans, excluding age factor
    Continuous.registerCallback( "alpha", ctsAlpha.str(), MakeDelegate( this, &VectorModel::ctsCbAlpha ) );
    Continuous.registerCallback( "P_B", ctsPB.str(), MakeDelegate( this, &VectorModel::ctsCbP_B ) );
    Continuous.registerCallback( "P_C*P_D", ctsPCD.str(), MakeDelegate( this, &VectorModel::ctsCbP_CD ) );
//     Continuous.registerCallback( "mean insecticide content",
//         "\tmean insecticide content",
//         MakeDelegate( this, &VectorModel::ctsNetInsecticideContent ) );
    // Mean IRS insecticide across whole population
//     Continuous.registerCallback( "IRS insecticide content",
//         "\tIRS insecticide content",
//         MakeDelegate( this, &VectorModel::ctsIRSInsecticideContent ) );
    // Mean values of relative attractiveness, pre- and post-prandial survival
    // IRS-induced factors of mosquitoes (i.e. only the portion of deterrent
    // and killing effects attributable to IRS).
//     Continuous.registerCallback( "IRS effects", ctsIRSEffects.str(),
//         MakeDelegate( this, &VectorModel::ctsIRSEffects ) );
    Continuous.registerCallback( "resource availability", ctsRA.str(),
        MakeDelegate( this, &VectorModel::ctsCbResAvailability ) );
    Continuous.registerCallback( "resource requirements", ctsRR.str(),
        MakeDelegate( this, &VectorModel::ctsCbResRequirements ) );
}
VectorModel::~VectorModel () {
}

void VectorModel::init2 () {
    double sumRelativeAvailability = 0.0;
    foreach(const Host::Human& human, sim::humanPop().crange()) {
        sumRelativeAvailability +=
                human.perHostTransmission.relativeAvailabilityAge (human.age(sim::now()).inYears());
    }
    int popSize = sim::humanPop().size();
    // value should be unimportant when no humans are available, though inf/nan is not acceptable
    double meanPopAvail = 1.0;
    if( popSize > 0 ){
        meanPopAvail = sumRelativeAvailability / popSize;    // mean-rel-avail
    }
    
    for(size_t i = 0; i < numSpecies; ++i) {
        species[i].init2 (i, meanPopAvail);
    }
    simulationMode = forcedEIR;   // now we should be ready to start
}

void VectorModel::initVectorInterv( const scnXml::Description::AnophelesSequence& list,
        size_t instance, const string& name )
{
    SpeciesIndexChecker checker( name, speciesIndex );
    foreach( const scnXml::VectorSpeciesIntervention& anoph, list ){
        const string& mosq = anoph.getMosquito();
        species[checker.getIndex(mosq)].initVectorInterv( anoph, instance );
    }
    checker.checkNoneMissed();
}
void VectorModel::initVectorTrap( const scnXml::VectorTrap::DescriptionSequence list,
        size_t instance, const scnXml::VectorTrap::NameOptional name_opt )
{
    stringstream name_ss;
    if( name_opt.present() ) name_ss << name_opt.get();
    else name_ss << "vector trap intervention " << (instance+1);
    string name = name_ss.str();
    SpeciesIndexChecker checker( name, speciesIndex );
    foreach( const scnXml::Description1& anoph, list ){
        const string& mosq = anoph.getMosquito();
        species[checker.getIndex(mosq)].initVectorTrap( anoph, instance );
    }
    checker.checkNoneMissed();
}


void VectorModel::scaleEIR (double factor) {
    for( size_t i = 0; i < numSpecies; ++i )
        species[i].scaleEIR( factor );
    vectors::scale( initialisationEIR, factor );
    annualEIR = vectors::sum( initialisationEIR );
}
#if 0
void VectorModel::scaleXML_EIR (scnXml::EntoData& ed, double factor) const {
    // XML values are exponentiated; so we add some factor to existing a0 values:
    double add_to_a0 = std::log( factor );

    assert( ed.getVector().present() );
    scnXml::Vector::AnophelesSequence& anophelesList = ed.getVector().get().getAnopheles();

    for( scnXml::Vector::AnophelesIterator it = anophelesList.begin();
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
#endif


SimTime VectorModel::minPreinitDuration () {
    if ( interventionMode == forcedEIR ) {
        return SimTime::zero();
    }
    // Data is summed over 5 years; add an extra 50 for stabilization.
    // 50 years seems a reasonable figure from a few tests
    return SimTime::fromYearsI( 55 );
}
SimTime VectorModel::expectedInitDuration (){
    return SimTime::oneYear();
}

SimTime VectorModel::initIterate () {
    if( interventionMode != dynamicEIR ) {
        // allow forcing equilibrium mode like with non-vector model
        return SimTime::zero(); // no initialization to do
    }
    if( initIterations < 0 ){
        assert( interventionMode = dynamicEIR );
        simulationMode = dynamicEIR;
        return SimTime::zero();
    }
    
    ++initIterations;
    
    bool needIterate = false;
    for(size_t i = 0; i < numSpecies; ++i) {
        needIterate = needIterate || species[i].initIterate ();
    }
    if( needIterate == false ){
        initIterations = -1;
    }
    if( initIterations > 10 ){
        throw TRACED_EXCEPTION("Transmission warmup exceeded 10 iterations!",util::Error::VectorWarmup);
    }
    
    // Time to let parameters settle after each iteration. I would expect one year
    // to be enough (but I may be wrong).
    if( needIterate )
        // stabilization + 5 years data-collection time:
        return SimTime::oneYear() + SimTime::fromYearsI(5);
    else
        return SimTime::oneYear();
}

double VectorModel::calculateEIR(Host::Human& human, double ageYears,
        vector<double>& EIR)
{
    PerHost& host = human.perHostTransmission;
    host.update( human );
    if (simulationMode == forcedEIR){
        double eir = initialisationEIR[sim::ts0().moduloYearSteps()] *
                host.relativeAvailabilityHetAge (ageYears);
        EIR.assign( 1, eir );
        return eir;
    }else{
        assert( simulationMode == dynamicEIR );
        EIR.assign( WithinHost::Genotypes::N(), 0.0 );
        for(size_t i = 0; i < numSpecies; ++i) {
            species[i].calculateEIR( i, host, EIR );
        }
        vectors::scale( EIR, host.relativeAvailabilityAge (ageYears) );
        return vectors::sum( EIR );
    }
}


// Every Global::interval days:
void VectorModel::vectorUpdate () {
    vector2D<double> popProbTransmission;
    popProbTransmission.resize( sim::humanPop().size(), WithinHost::Genotypes::N() );
    size_t i = 0;
    foreach(const Host::Human& human, sim::humanPop().crange()) {
        const double tbvFac = human.getVaccine().getFactor( interventions::Vaccine::TBV );
        WithinHost::WHInterface& whm = *human.withinHostModel;
        double sumX;
        const double pTrans = whm.probTransmissionToMosquito( tbvFac, &sumX );
        if( WithinHost::Genotypes::N() == 1 ) popProbTransmission.at(i,0) = pTrans;
        else for( size_t g = 0; g < WithinHost::Genotypes::N(); ++g ){
            const double k = whm.probTransGenotype( pTrans, sumX, g );
            assert( (boost::math::isfinite)(k) );
            popProbTransmission.at(i,g) = k;
        }
        i += 1;
    }
    for(size_t i = 0; i < numSpecies; ++i){
        species[i].advancePeriod (popProbTransmission, i, simulationMode == dynamicEIR);
    }
}
void VectorModel::update() {
    TransmissionModel::updateKappa();
}

const string vec_mode_err = "vector interventions can only be used in "
            "dynamic transmission mode (mode=\"dynamic\")";
const map<string,size_t>& VectorModel::getSpeciesIndexMap(){
    if( interventionMode != dynamicEIR ){
        throw xml_scenario_error(vec_mode_err);
    }
    return speciesIndex;
}

void VectorModel::deployVectorPopInterv (size_t instance) {
    if( interventionMode != dynamicEIR ){
        throw xml_scenario_error(vec_mode_err);
    }
    for( vector<AnophelesModel>::iterator it = species.begin(); it != species.end(); ++it ){
        it->deployVectorPopInterv(instance);
    }
}
void VectorModel::deployVectorTrap( size_t instance, double number, SimTime lifespan ){
    if( interventionMode != dynamicEIR ){
        throw xml_scenario_error(vec_mode_err);
    }
    for( vector<AnophelesModel>::iterator it = species.begin(); it != species.end(); ++it ){
        it->deployVectorTrap( instance, number, lifespan );
    }
}
void VectorModel::uninfectVectors() {
    for(size_t i = 0; i < numSpecies; ++i)
        species[i].uninfectVectors();
}

void VectorModel::summarize () {
    TransmissionModel::summarize ();
    
    for(size_t i = 0; i < numSpecies; ++i){
        species[i].summarize( i );
    }
}


void VectorModel::checkpoint (istream& stream) {
    TransmissionModel::checkpoint (stream);
    initIterations & stream;
    species & stream;
}
void VectorModel::checkpoint (ostream& stream) {
    TransmissionModel::checkpoint (stream);
    initIterations & stream;
    species & stream;
}

}
}
