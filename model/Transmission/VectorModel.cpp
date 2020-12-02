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
#include "mon/Continuous.h"
#include "util/vectors.h"
#include "util/ModelOptions.h"
#include "util/SpeciesIndexChecker.h"
#include "util/StreamValidator.h"
#include "Transmission/Anopheles/SimpleMPDAnophelesModel.h"

#include <fstream>
#include <map>
#include <cmath>
#include <set>

namespace OM {
namespace Transmission {
using namespace OM::util;
using Anopheles::PerHostAnophParams;
using Anopheles::AnophelesModel;
using Anopheles::SimpleMPDAnophelesModel;

/** A map of anopheles species/variant name to an index in species.
 *
 * When the main ento data is read from XML, each anopheles section is read
 * into one index of the species array. The name and this index are added
 * here.
 * 
 * Other data read from XML should look up the name here and use the index
 * found. Doesn't need checkpointing. */
map<string,size_t> speciesIndex;


void VectorModel::ctsCbN_v0 (ostream& stream) {
    for(size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i]->getLastN_v0();
}
void VectorModel::ctsCbP_A (ostream& stream) {
    for(size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i]->getLastVecStat(Anopheles::PA);
}
void VectorModel::ctsCbP_Amu (ostream& stream) {
    for(size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i]->getLastVecStat(Anopheles::PAmu);
}
void VectorModel::ctsCbP_A1 (ostream& stream) {
    for(size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i]->getLastVecStat(Anopheles::PA1);
}
void VectorModel::ctsCbP_Ah (ostream& stream) {
    for(size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i]->getLastVecStat(Anopheles::PAh);
}
void VectorModel::ctsCbP_Amu (ostream& stream) {
    for(size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i].getLastVecStat(Anopheles::PAmu);
}
void VectorModel::ctsCbP_A1 (ostream& stream) {
    for(size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i].getLastVecStat(Anopheles::PA1);
}
void VectorModel::ctsCbP_Ah (ostream& stream) {
    for(size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i].getLastVecStat(Anopheles::PAh);
}
void VectorModel::ctsCbP_df (ostream& stream) {
    for(size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i]->getLastVecStat(Anopheles::PDF);
}
void VectorModel::ctsCbP_dif (ostream& stream) {
    for(size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i]->getLastVecStat(Anopheles::PDIF);
}
void VectorModel::ctsCbN_v (ostream& stream) {
    for(size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i]->getLastVecStat(Anopheles::NV);
}
void VectorModel::ctsCbO_v (ostream& stream) {
    for(size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i]->getLastVecStat(Anopheles::OV);
}
void VectorModel::ctsCbS_v (ostream& stream) {
    for(size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i]->getLastVecStat(Anopheles::SV);
}
void VectorModel::ctsCbAlpha (const Population& population, ostream& stream){
    for( size_t i = 0; i < speciesIndex.size(); ++i){
        double total = 0.0;
        for(Population::ConstIter iter = population.cbegin(); iter != population.cend(); ++iter) {
            total += iter->perHostTransmission.entoAvailabilityFull( i, iter->age(sim::now()).inYears() );
        }
        stream << '\t' << total / population.size();
    }
}
void VectorModel::ctsCbP_B (const Population& population, ostream& stream){
    for( size_t i = 0; i < speciesIndex.size(); ++i){
        double total = 0.0;
        for(Population::ConstIter iter = population.cbegin(); iter != population.cend(); ++iter) {
            total += iter->perHostTransmission.probMosqBiting( i );
        }
        stream << '\t' << total / population.size();
    }
}
void VectorModel::ctsCbP_CD (const Population& population, ostream& stream){
    for( size_t i = 0; i < speciesIndex.size(); ++i){
        double total = 0.0;
        for(Population::ConstIter iter = population.cbegin(); iter != population.cend(); ++iter) {
            total += iter->perHostTransmission.probMosqResting( i );
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
//     for( size_t i = 0; i < speciesIndex.size(); ++i ){
//         const interventions::IRSAnophelesParams& params = species[i]->getHumanBaseParams().irs;
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
    for(size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i]->getResAvailability();
}
void VectorModel::ctsCbResRequirements (ostream& stream) {
    for(size_t i = 0; i < speciesIndex.size(); ++i)
        stream << '\t' << species[i]->getResRequirements();
}

const string& reverseLookup (const map<string,size_t>& m, size_t i) {
    for( auto it = m.begin(); it != m.end(); ++it ) {
        if ( it->second == i )
            return it->first;
    }
    throw TRACED_EXCEPTION_DEFAULT( "reverseLookup: key not found" );        // shouldn't ever happen
}

bool anophelesCompare(const scnXml::AnophelesParams& a1, const scnXml::AnophelesParams& a2)
{
    return (a1.getSeasonality().getAnnualEIR().get()>a2.getSeasonality().getAnnualEIR().get()); 
}

VectorModel::VectorModel (
                          const scnXml::Entomology& entoData,
                          const scnXml::Vector vectorData, int populationSize) :
    TransmissionModel( entoData, WithinHost::Genotypes::N() ),
    m_rng(util::master_RNG), initIterations(0)
{
    // Each item in the AnophelesSequence represents an anopheles species.
    // TransmissionModel::createTransmissionModel checks length of list >= 1.
    scnXml::Vector::AnophelesSequence anophelesList = vectorData.getAnopheles();
    const scnXml::Vector::NonHumanHostsSequence nonHumansList = vectorData.getNonHumanHosts();

//     map<string, double> nonHumanHostPopulations;
//     for(size_t i = 0; i<nonHumansList.size(); i++) {
//         nonHumanHostPopulations[nonHumansList[i].getName()] = nonHumansList[i].getNumber();
//     }

    size_t numSpecies = anophelesList.size();
    if (numSpecies < 1)
        throw util::xml_scenario_error ("Can't use Vector model without data for at least one anopheles species!");
    
    sort(anophelesList.begin(), anophelesList.end(), anophelesCompare);

    PerHostAnophParams::initReserve (numSpecies);
    species.resize (numSpecies);

//     for(size_t i = 0; i < numSpecies; ++i) {
//         auto elt = anophelesList[i];
//         PerHostAnophParams::init(elt.getMosq());
//         string name = species[i]->initialise (i, elt,
//                                              initialisationEIR,
// //                                              nonHumanHostPopulations,
//                                              populationSize);
//         speciesIndex[name] = i;
//     }

    for(size_t i = 0; i < numSpecies; ++i)
    {
        auto elt = anophelesList[i];
 
        std::unique_ptr<Anopheles::AnophelesModel> anophModel;

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
            const scnXml::AnophelesParams::SimpleMPDOptional& simpleMPDOpt = elt.getSimpleMPD();

            if (!simpleMPDOpt.present())
                throw util::xml_scenario_error("VECTOR_SIMPLE_MPD_MODEL: requires <simpleMPD> element with "
                                               "model parameters for each anopheles species");
            anophModel = std::unique_ptr<Anopheles::SimpleMPDAnophelesModel>(new Anopheles::SimpleMPDAnophelesModel(simpleMPDOpt.get()));
            // anophModel = std::unique_ptr<Anopheles::AnophelesModel>(new Anopheles::AnophelesModel());        
        }
        else
            anophModel = std::unique_ptr<Anopheles::AnophelesModel>(new Anopheles::AnophelesModel());

        PerHostAnophParams::init(elt.getMosq());

        species[i] = std::move(anophModel);
        species[i]->initialise(i, elt, initialisationEIR, populationSize);

        string name = elt.getMosquito();
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
    ostringstream ctsNv0, ctsPA, ctsPAmu, ctsPA1, ctsPAh,
        ctsPdf, ctsPdif,
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
        ctsPAmu<<"\tP_Amu("<<name<<")";
        ctsPA1<<"\tP_A1("<<name<<")";
        ctsPAh<<"\tP_Ah("<<name<<")";
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
    using mon::Continuous;
    Continuous.registerCallback( "N_v0", ctsNv0.str(), MakeDelegate( this, &VectorModel::ctsCbN_v0 ) );
    Continuous.registerCallback( "P_A", ctsPA.str(), MakeDelegate( this, &VectorModel::ctsCbP_A ) );
    Continuous.registerCallback( "P_Amu", ctsPAmu.str(), MakeDelegate( this, &VectorModel::ctsCbP_Amu ) );
    Continuous.registerCallback( "P_A1", ctsPA1.str(), MakeDelegate( this, &VectorModel::ctsCbP_A1 ) );
    Continuous.registerCallback( "P_Ah", ctsPAh.str(), MakeDelegate( this, &VectorModel::ctsCbP_Ah ) );
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

void VectorModel::init2 (const Population& population) {
    SimTime data_save_len = SimTime::oneDay();  // we don't need to save anything at first
    saved_sum_avail.assign( data_save_len, speciesIndex.size(), 0.0 );
    saved_sigma_df.assign( data_save_len, speciesIndex.size(), 0.0 );
    saved_sigma_dif.assign( data_save_len, speciesIndex.size(), WithinHost::Genotypes::N(), 0.0 );
    saved_sigma_dff.assign( speciesIndex.size(), 0.0 );
    
    double sumRelativeAvailability = 0.0;
    for(const Host::Human& human : population.getHumans()) {
        sumRelativeAvailability +=
                human.perHostTransmission.relativeAvailabilityAge (human.age(sim::now()).inYears());
    }
    int popSize = population.size();
    // value should be unimportant when no humans are available, though inf/nan is not acceptable
    double meanPopAvail = 1.0;
    if( popSize > 0 ){
        meanPopAvail = sumRelativeAvailability / popSize;    // mean-rel-avail
    }
    
    for(size_t i = 0; i < speciesIndex.size(); ++i) {
        double sum_avail = 0.0;
        double sigma_f = 0.0;
        double sigma_df = 0.0;
        double sigma_dff = 0.0;
        
        for(const Host::Human& human : population.getHumans()) {
            const OM::Transmission::PerHost& host = human.perHostTransmission;
            double prod = host.entoAvailabilityFull (i, human.age(sim::now()).inYears());
            sum_avail += prod;
            prod *= host.probMosqBiting(i);
            sigma_f += prod;
            sigma_df += prod * host.probMosqResting(i);
            sigma_dff += prod * host.probMosqResting(i) * host.relMosqFecundity(i);
        }
        
        species[i]->init2 (population.size(), meanPopAvail, sum_avail, sigma_f, sigma_df, sigma_dff);
    }
    simulationMode = forcedEIR;   // now we should be ready to start
}

void VectorModel::initVectorInterv( const scnXml::Description::AnophelesSequence& list,
        size_t instance, const string& name )
{
    SpeciesIndexChecker checker( name, speciesIndex );
    for( const scnXml::VectorSpeciesIntervention& anoph : list ){
        const string& mosq = anoph.getMosquito();
        species[checker.getIndex(mosq)]->initVectorInterv( anoph, instance );
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
    for( const scnXml::Description1& anoph : list ){
        const string& mosq = anoph.getMosquito();
        species[checker.getIndex(mosq)]->initVectorTrap( anoph, instance );
    }
    checker.checkNoneMissed();
}
void VectorModel::initNonHumanHostsInterv( const scnXml::Description2::AnophelesSequence list,
        const scnXml::DecayFunction& decay, size_t instance, const string& name )
{
    SpeciesIndexChecker checker( name, speciesIndex );
    for( const scnXml::NonHumanHostsSpeciesIntervention& anoph : list ){
        const string& mosq = anoph.getMosquito();
        species[checker.getIndex(mosq)]->initNonHumanHostsInterv( anoph, decay, instance, name );
    }
    checker.checkNoneMissed();
}
void VectorModel::initAddNonHumanHostsInterv( const scnXml::Description3::AnophelesSequence list, const string& name )
{
    SpeciesIndexChecker checker( name, speciesIndex );
    for( const scnXml::NonHumanHostsVectorSpecies& anoph : list ){
        const string& mosq = anoph.getMosquito();
        species[checker.getIndex(mosq)]->initAddNonHumanHostsInterv( anoph, name );
    }
    checker.checkNoneMissed();
}

void VectorModel::scaleEIR (double factor) {
    for( size_t i = 0; i < speciesIndex.size(); ++i )
        species[i]->scaleEIR( factor );
    vectors::scale( initialisationEIR, factor );
    annualEIR = vectors::sum( initialisationEIR );
}
#if 0
void VectorModel::scaleXML_EIR (scnXml::EntoData& ed, double factor) const {
    // XML values are exponentiated; so we add some factor to existing a0 values:
    double add_to_a0 = std::log( factor );

    assert( ed.getVector().present() );
    scnXml::Vector::AnophelesSequence& anophelesList = ed.getVector().get().getAnopheles();

    for( auto it = anophelesList.begin();
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
    
    // This function is called repeatedly until vector initialisation is
    // complete (signalled by returning 0).
    int initState = 0;
    if( initIterations == 0 && saved_sum_avail.size1() == SimTime::oneDay() ) {
        // First time called: we need to do some data collection
        initState = 1;
    } else if( initIterations >= 0 ){
        // Next, while generated EIR is not close to that required,
        // try adjusting emergence parameters, do a stabilisation phase then data collection phase, then repeat.
        initState = 2;
    } else if( initIterations < 0 ){
        // Finally, we're done.
        initState = 3;
    }
    
    if( initState == 1 || initState == 3 ){
        // When starting the iteration phase, we switch to five years worth of data; otherwise we only keep one day.
        SimTime data_save_len = /*initState == 1 ? SimTime::fromYearsI(5) :*/ SimTime::oneDay();
        saved_sum_avail.assign( data_save_len, speciesIndex.size(), 0.0 );
        saved_sigma_df.assign( data_save_len, speciesIndex.size(), 0.0 );
        saved_sigma_dif.assign( data_save_len, speciesIndex.size(), WithinHost::Genotypes::N(), 0.0 );
        
        if( initState == 3 ){
            //TODO: we should perhaps check that EIR gets reproduced correctly?
            simulationMode = dynamicEIR;
            return SimTime::zero();
        }
    }
    
    ++initIterations;
    if( initIterations > 30 ){
        throw TRACED_EXCEPTION("Transmission warmup exceeded 30 iterations!",util::Error::VectorWarmup);
    }
    
    bool needIterate = false;
    for(size_t i = 0; i < speciesIndex.size(); ++i) {
        needIterate = needIterate || species[i]->initIterate ();
    }
    
    if( needIterate ){
        // stabilization + 5 years data-collection time:
        return SimTime::oneYear() + SimTime::fromYearsI(5);
    } else {
        // One year stabilisation, then we're finished:
        initIterations = -1;
        return SimTime::oneYear();
    }
}

void VectorModel::calculateEIR(Host::Human& human, double ageYears,
        vector<double>& EIR) const
{
    auto ag = human.monAgeGroup().i();
    auto cs = human.cohortSet();
    PerHost& host = human.perHostTransmission;
    host.update( human );
    if (simulationMode == forcedEIR){
        double eir = initialisationEIR[sim::ts0().moduloYearSteps()] *
                host.relativeAvailabilityHetAge (ageYears);
        mon::reportStatMACGF( mon::MVF_INOCS, ag, cs, 0, eir );
        EIR.assign( 1, eir );
    }else{
        assert( simulationMode == dynamicEIR );
        EIR.assign( WithinHost::Genotypes::N(), 0.0 );
        const double ageFactor = host.relativeAvailabilityAge (ageYears);
        for(size_t i = 0; i < speciesIndex.size(); ++i) {
            const vector<double>& partialEIR = species[i]->getPartialEIR();
            
            assert( EIR.size() == partialEIR.size() );
            if ( (std::isnan)(vectors::sum(partialEIR)) ) {
                cerr<<"partialEIR is not a number; "<<i<<endl;
            }
            
            /* Calculates EIR per individual (hence N_i == 1).
             *
             * See comment in AnophelesModel::advancePeriod for method. */
            double entoFactor = ageFactor * host.availBite(i);
            for( size_t g = 0; g < EIR.size(); ++g ){
                auto eir = partialEIR[g] * entoFactor;
                mon::reportStatMACSGF( mon::MVF_INOCS, ag, cs, i, g, eir );
                EIR[g] += eir;
            }
        }
    }
}


// Every Global::interval days:
void VectorModel::vectorUpdate (const Population& population) {
    const size_t nGenotypes = WithinHost::Genotypes::N();
    SimTime popDataInd = mod_nn(sim::ts0(), saved_sum_avail.size1());
    vector<double> probTransmission;
    saved_sum_avail.assign_at1(popDataInd, 0.0);
    saved_sigma_df.assign_at1(popDataInd, 0.0);
    saved_sigma_dif.assign_at1(popDataInd, 0.0);
    saved_sigma_dff.assign( saved_sigma_dff.size(), 0.0 );
    
    for(const Host::Human& human : population.getHumans()) {
        const OM::Transmission::PerHost& host = human.perHostTransmission;
        WithinHost::WHInterface& whm = *human.withinHostModel;
        const double tbvFac = human.getVaccine().getFactor( interventions::Vaccine::TBV );
        
        probTransmission.assign( nGenotypes, 0.0 );
        double sumX = numeric_limits<double>::quiet_NaN();
        const double pTrans = whm.probTransmissionToMosquito( tbvFac, &sumX );
        if( nGenotypes == 1 ) probTransmission[0] = pTrans;
        else for( size_t g = 0; g < nGenotypes; ++g ){
            const double k = whm.probTransGenotype( pTrans, sumX, g );
            assert( (std::isfinite)(k) );
            probTransmission[g] = k;
        }
        
        for(size_t s = 0; s < speciesIndex.size(); ++s){
            //NOTE: calculate availability relative to age at end of time step;
            // not my preference but consistent with TransmissionModel::getEIR().
            //TODO: even stranger since probTransmission comes from the previous time step
            const double avail = host.entoAvailabilityFull (s, human.age(sim::ts1()).inYears());
            saved_sum_avail.at(popDataInd, s) += avail;
            const double df = avail
                    * host.probMosqBiting(s)
                    * host.probMosqResting(s);
            saved_sigma_df.at(popDataInd, s) += df;
            for( size_t g = 0; g < nGenotypes; ++g ){
                saved_sigma_dif.at(popDataInd, s, g) += df * probTransmission[g];
            }
            saved_sigma_dff[s] += df * host.relMosqFecundity(s);
        }
    }
    
    for(size_t s = 0; s < speciesIndex.size(); ++s){
        // Copy slice to new array:
        auto range = saved_sigma_dif.range_at12(popDataInd, s);
        sigma_dif_species.assign(range.first, range.second);
        
        species[s]->advancePeriod (saved_sum_avail.at(popDataInd, s),
                saved_sigma_df.at(popDataInd, s),
                sigma_dif_species,
                saved_sigma_dff[s],
                simulationMode == dynamicEIR);
    }
}
void VectorModel::update(const Population& population) {
    TransmissionModel::updateKappa(population);
}

const string vec_mode_err = "vector interventions can only be used in "
            "dynamic transmission mode (mode=\"dynamic\")";
const map<string,size_t>& VectorModel::getSpeciesIndexMap() {
    return speciesIndex;
}

void VectorModel::deployVectorPopInterv (size_t instance) {
    if( interventionMode != dynamicEIR ){
        throw xml_scenario_error(vec_mode_err);
    }
    for( auto it = species.begin(); it != species.end(); ++it ){
        (*it)->deployVectorPopInterv(m_rng, instance);
    }
}
void VectorModel::deployVectorTrap( size_t instance, double number, SimTime lifespan ){
    if( interventionMode != dynamicEIR ){
        throw xml_scenario_error(vec_mode_err);
    }
    for(size_t i = 0; i < speciesIndex.size(); ++i) {
        species[i]->deployVectorTrap( m_rng, i, instance, number, lifespan );
    }
}
void VectorModel::deployNonHumanHostsInterv( size_t instance, string name ){
    if( interventionMode != dynamicEIR ){
        throw xml_scenario_error(vec_mode_err);
    }
    for(size_t i = 0; i < speciesIndex.size(); ++i) {
        species[i]->deployNonHumanHostsInterv( m_rng, i, instance, name);
    }
}
void VectorModel::deployAddNonHumanHosts(string name, double popSize, SimTime lifespan)
{
    if( interventionMode != dynamicEIR ){
        throw xml_scenario_error(vec_mode_err);
    }
    for(size_t i = 0; i < speciesIndex.size(); ++i) {
        species[i]->deployAddNonHumanHosts( m_rng, i, name, popSize, lifespan);
    }
}
void VectorModel::deployNonHumanHostsInterv( size_t instance, string name ){
    if( interventionMode != dynamicEIR ){
        throw xml_scenario_error(vec_mode_err);
    }
    for(size_t i = 0; i < speciesIndex.size(); ++i) {
        species[i].deployNonHumanHostsInterv( m_rng, i, instance, name);
    }
}
void VectorModel::deployAddNonHumanHosts(string name, double popSize, SimTime lifespan)
{
    if( interventionMode != dynamicEIR ){
        throw xml_scenario_error(vec_mode_err);
    }
    for(size_t i = 0; i < speciesIndex.size(); ++i) {
        species[i].deployAddNonHumanHosts( m_rng, i, name, popSize, lifespan);
    }
}
void VectorModel::uninfectVectors() {
    for(size_t i = 0; i < speciesIndex.size(); ++i)
        species[i]->uninfectVectors();
}

void VectorModel::summarize () {
    TransmissionModel::summarize ();
    
    for(size_t i = 0; i < speciesIndex.size(); ++i){
        species[i]->summarize( i );
    }
}


void VectorModel::checkpoint (istream& stream) {
    TransmissionModel::checkpoint (stream);
    m_rng.checkpoint(stream);
    initIterations & stream;
    // species & stream;
    for(auto &s : species)
        s->checkpoint(stream);
    saved_sum_avail & stream;
    saved_sigma_df & stream;
    saved_sigma_dif & stream;
    saved_sigma_dff & stream;
}
void VectorModel::checkpoint (ostream& stream) {
    TransmissionModel::checkpoint (stream);
    m_rng.checkpoint(stream);
    initIterations & stream;
    for(auto &s : species)
        s->checkpoint(stream);
    // species & stream;
    saved_sum_avail & stream;
    saved_sigma_df & stream;
    saved_sigma_dif & stream;
    saved_sigma_dff & stream;
}

}
}
