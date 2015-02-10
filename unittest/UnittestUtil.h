/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2014 Swiss Tropical and Public Health Institute
 Copyright (C) 2005-2014 Liverpool School Of Tropical Medicine
 
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

// Utility for unittests, which is granted "friend" access where necessary in the model.

#ifndef Hmod_UnittestUtil
#define Hmod_UnittestUtil

#include "Global.h"
#include "util/ModelOptions.h"

#include "Host/Human.h"
#include "PkPd/LSTMModel.h"
#include "PkPd/Drug/LSTMDrugType.h"
#include "PkPd/LSTMTreatments.h"
#include "WithinHost/Infection/Infection.h"
#include "WithinHost/WHFalciparum.h"
#include "WithinHost/Infection/MolineauxInfection.h"
#include "WithinHost/Genotypes.h"
#include "mon/management.h"

#include "schema/scenario.h"

using namespace OM;
using namespace WithinHost;

using ::OM::util::ModelOptions;

namespace OM {
    namespace WithinHost {
        extern bool opt_common_whm;
    }
}

namespace dummyXML{
    scnXml::DemogAgeGroup demAgeGroup(
        numeric_limits<double>::quiet_NaN() /* lower bound */ );
    scnXml::Demography demography( 
        demAgeGroup,
        "dummy" /* name */,
        0 /* pop size */,
        90.0 /*max human age*/ );
    
    scnXml::MonitoringOptions survOpts;
    scnXml::Surveys surveys;
    scnXml::MonAgeGroup monAgeGroup(
        0.0 /* lower bound */ );
    scnXml::Monitoring monitoring(
        survOpts,
        surveys,
        monAgeGroup,
        "dummy" /* name */ );
    
    scnXml::Interventions interventions( "dummy" /* name */ );
    
    scnXml::AgeGroupValues cfr, pSeq;
    scnXml::HealthSystem healthSystem( cfr, pSeq );
    
    scnXml::Entomology entomology( "dummy" /* name */,
        "dummy" /* mode */ );
    
    scnXml::OptionSet modelOpts;
    scnXml::Clinical modelClinical( "dummy" /* HS memory */ );
    scnXml::AgeGroupValues modelHumanAvailMosq;
    scnXml::Human modelHuman( modelHumanAvailMosq );
    scnXml::Parameters modelParams(
        0 /* interval */,
        0 /* iseed */,
        "dummy" /* latentP */ );
    scnXml::Model model( modelOpts, modelClinical, modelHuman, modelParams );
    
    scnXml::Scenario scenario(
        demography,
        monitoring,
        interventions,
        healthSystem,
        entomology,
        model,
        0 /* schema version */,
        "dummy" /* name */ );
}

class UnittestUtil {
public:
    static void initTime(int daysPerStep){
        dummyXML::modelParams.setInterval( daysPerStep );
        dummyXML::model.setParameters(dummyXML::modelParams);
        dummyXML::scenario.setModel(dummyXML::model);
        sim::init( dummyXML::scenario );
        
        // we could just use zero, but we may spot more errors by using some weird number
        sim::time0 = sim::fromYearsN(83.2591);
        sim::time1 = sim::time0;
#ifndef NDEBUG
        sim::in_update = true;  // may not always be correct but we're more interested in getting around this check than using it in unit tests
#endif
    }
    static void incrTime(SimTime incr){
        //NOTE: for unit tests, we do not differentiate between time0 and time1
        sim::time0 += incr;
        sim::time1 = sim::time0;
    }
    
    static const scnXml::Parameters& prepareParameters(){
        if( dummyXML::modelParams.getParameter().size() == 0 ){
            dummyXML::modelParams.getParameter().push_back( scnXml::Parameter( 15, 0.177378570987455 ) );
            dummyXML::modelParams.getParameter().push_back( scnXml::Parameter( 34, 4.7601 ) );
            dummyXML::modelParams.getParameter().push_back( scnXml::Parameter( 35, 0.5008 ) );
            dummyXML::modelParams.getParameter().push_back( scnXml::Parameter( 36, 2.2736 ) );
            dummyXML::modelParams.getParameter().push_back( scnXml::Parameter( 37, 0.2315 ) );
        }
        return dummyXML::modelParams;
    }
    
    // Initialise surveys, to the minimum required not to crash
    static void initSurveys(){
        diagnostics::clear();
        Parameters parameters( prepareParameters() );
        dummyXML::surveys.setDetectionLimit( numeric_limits<double>::quiet_NaN() );
        dummyXML::monitoring.setSurveys( dummyXML::surveys );
        mon::initSurveyTimes( parameters, dummyXML::scenario, dummyXML::monitoring );
    }
    // Parameterise standard diagnostics
    static void setDiagnostics(){
        // note that this is only ever called after initSurveys(), thus we
        // don't need to call diagnostics::clear() (and shouldn't because it
        // would leave a dangling pointer in Survey::m_diagnostic)
        scnXml::Diagnostic microscopy( "microscopy" );
        scnXml::Stochastic microscopyParams( 20, 0.75 );
        microscopy.setStochastic( microscopyParams );
        scnXml::Diagnostic rdt( "RDT" );
        scnXml::Stochastic rdtParams( 50, 0.942 );
        rdt.setStochastic( rdtParams );
        scnXml::Diagnostics diagsElt;
        diagsElt.getDiagnostic().push_back( microscopy );
        diagsElt.getDiagnostic().push_back( rdt );
        Parameters parameters( prepareParameters() );
        dummyXML::scenario.setDiagnostics(diagsElt);
        dummyXML::surveys.setDetectionLimit(numeric_limits<double>::quiet_NaN());
        dummyXML::monitoring.setSurveys(dummyXML::surveys);
        dummyXML::scenario.setMonitoring(dummyXML::monitoring);
        diagnostics::init( parameters, dummyXML::scenario );
    }
    
    static void PkPdSuiteSetup () {
        ModelOptions::reset();
        WithinHost::Genotypes::initSingle();

        //Note: we fudge this call since it's not so easy to falsely initialize         scenario element.
        //PkPdModel::init ();

        // Drugs:
        scnXml::Drugs drugs;
        // Artemether
        scnXml::Phenotype phenotypeAR ( 27.6 /* max_killing_rate */, 0.0023 /* IC50 */, 
            4.0 /* slope */ );
        scnXml::PD pdAR;
        pdAR.getPhenotype().push_back (phenotypeAR);
        scnXml::PK pkAR ( 1e-17 /* negligible_concentration */, 0.1750372 /* half_life */, 
            17.4 /* vol_dist */ );
        scnXml::PKPDDrug drugAR ( pdAR, pkAR, "AR" /* abbrev */ );
        
        drugs.getDrug().push_back (drugAR);
        
        // Artesunate
        scnXml::Phenotype phenotypeAS ( 27.6 /* max_killing_rate */, 0.0016 /* IC50 */, 
            4.0 /* slope */ );
        scnXml::PD pdAS;
        pdAS.getPhenotype().push_back (phenotypeAS);
        scnXml::PK pkAS ( 1e-17 /* negligible_concentration */, 0.04175585 /* half_life */, 
            2.75 /* vol_dist */ );
        scnXml::PKPDDrug drugAS ( pdAS, pkAS, "AS" /* abbrev */ );
        
        drugs.getDrug().push_back (drugAS);
                
        // Chloroquine
        scnXml::Phenotype phenotypeCQ ( 3.45 /* max_killing_rate */, 0.02 /* IC50 */, 
            1.6 /* slope */ );
        scnXml::PD pdCQ;
        pdCQ.getPhenotype().push_back (phenotypeCQ);
        scnXml::PK pkCQ ( 0.00036 /* negligible_concentration */, 30.006 /* half_life */, 
            300 /* vol_dist */ );
        scnXml::PKPDDrug drugCQ ( pdCQ, pkCQ, "CQ" /* abbrev */ );
        
        drugs.getDrug().push_back (drugCQ);
        
        // Dihydroartemisinin
        scnXml::Phenotype phenotypeDHA ( 27.6 /* max_killing_rate */, 0.009 /* IC50 */, 
            4.0 /* slope */ );
        scnXml::PD pdDHA;
        pdDHA.getPhenotype().push_back (phenotypeDHA);
        scnXml::PK pkDHA ( 1e-17 /* negligible_concentration */, 0.03500743 /* half_life */, 
            1.49 /* vol_dist */ );
        scnXml::PKPDDrug drugDHA ( pdDHA, pkDHA, "DHA" /* abbrev */ );
        
        drugs.getDrug().push_back (drugDHA);
                
        // Lumefantrine
        scnXml::Phenotype phenotypeLF ( 3.45 /* max_killing_rate */, 0.032 /* IC50 */, 
            4.0 /* slope */ );
        scnXml::PD pdLF;
        pdLF.getPhenotype().push_back (phenotypeLF);
        scnXml::PK pkLF ( 0.00032 /* negligible_concentration */, 4.332 /* half_life */, 
            21 /* vol_dist */ );
        scnXml::PKPDDrug drugLF ( pdLF, pkLF, "LF" /* abbrev */ );
        
        drugs.getDrug().push_back (drugLF);
        
        // Mefloquine
        scnXml::Phenotype phenotypeMQ ( 3.45 /* max_killing_rate */, 0.027 /* IC50 */, 
            5.0 /* slope */ );
        scnXml::PD pdMQ;
        pdMQ.getPhenotype().push_back (phenotypeMQ);
        scnXml::PK pkMQ ( 0.005 /* negligible_concentration */, 13.078 /* half_life */, 
            20.8 /* vol_dist */ );
        scnXml::PKPDDrug drugMQ ( pdMQ, pkMQ, "MQ" /* abbrev */ );
        drugs.getDrug().push_back (drugMQ);
        
        // Piperaquine
        scnXml::Phenotype phenotypePPQ ( 3.45 /* max_killing_rate */, 0.088 /* IC50 */, 
            6.0 /* slope */ );
        scnXml::PD pdPPQ;
        pdPPQ.getPhenotype().push_back (phenotypePPQ);
        scnXml::PK pkPPQ ( 0.005 /* negligible_concentration */, 23.105 /* half_life */, 
            150 /* vol_dist */ );
        scnXml::PKPDDrug drugPPQ ( pdPPQ, pkPPQ, "PPQ" /* abbrev */ );
        
        drugs.getDrug().push_back (drugPPQ);
        
        
        PkPd::LSTMDrugType::init (drugs);
        
        // Treatments
        scnXml::PKPDSchedule sched1("sched1");
        sched1.getMedicate().push_back(
            scnXml::PKPDMedication("MQ", 6 /*mg*/,0 /*hour*/ )
        );
        
        scnXml::PKPDSchedule sched2("sched2");
        sched2.getMedicate().push_back(
            scnXml::PKPDMedication("MQ", 2 /*mg*/, 0 /*hour*/));
        sched2.getMedicate().push_back(
            scnXml::PKPDMedication("MQ", 5 /*mg*/, 12 /*hour*/));
        
        // a very basic dosage table, so that we can test it does what's expected
        scnXml::PKPDDosages dosage1("dosage1");
        dosage1.getAge().push_back(scnXml::PKPDDosageRange(0 /*age lb*/,1 /*mult*/));
        dosage1.getAge().push_back(scnXml::PKPDDosageRange(5 /*age lb*/,5 /*mult*/));
        
        scnXml::Treatments treatments;
        treatments.getSchedule().push_back(sched1);
        treatments.getSchedule().push_back(sched2);
        treatments.getDosages().push_back(dosage1);
        PkPd::LSTMTreatments::init(treatments);
    }
    
    // For when infection parameters shouldn't be used; enforce by setting to NaNs.
    // But do set latentP.
    static void Infection_init_latentP_and_NaN () {
	Infection::latentP = sim::fromDays(15);
	Infection::invCumulativeYstar = numeric_limits<double>::quiet_NaN();
	Infection::invCumulativeHstar = numeric_limits<double>::quiet_NaN();
	Infection::alpha_m = numeric_limits<double>::quiet_NaN();
	Infection::decayM = numeric_limits<double>::quiet_NaN();
    }
    static void Infection_init_5day () {
	// Note: these values were pulled from one source and shouldn't be taken as authoritative
	Infection::latentP = sim::fromDays(15);
	Infection::invCumulativeYstar = 1.0 / 68564384.7102;
	Infection::invCumulativeHstar = 1.0 / 71.676733;
	Infection::alpha_m = 1.0 - exp(- 2.411434);
	Infection::decayM = 2.717773;
    }
    
    static void DescriptiveInfection_init () {
        ModelOptions::reset();
    }
    
    static void EmpiricalWHM_setup () {
        ModelOptions::reset();
        ModelOptions::set(util::EMPIRICAL_WITHIN_HOST_MODEL);
        OM::WithinHost::opt_common_whm = true;
    }
    
    static void MolineauxWHM_setup( const std::string& mode, bool repl_gamma ){
        ModelOptions::reset();
        ModelOptions::set(util::MOLINEAUX_WITHIN_HOST_MODEL);
        OM::WithinHost::opt_common_whm = true;
        if( mode == "original" ){
        }else if( mode == "1st_max_gamma" ){
            ModelOptions::set(util::FIRST_LOCAL_MAXIMUM_GAMMA);
        }else if( mode == "mean_dur_gamma" ){
            ModelOptions::set(util::MEAN_DURATION_GAMMA);
        }else if( mode == "1st_max_and_mean_dur_gamma" ){
            ModelOptions::set(util::FIRST_LOCAL_MAXIMUM_GAMMA);
            ModelOptions::set(util::MEAN_DURATION_GAMMA);
        }else if( mode == "pairwise" ){
            ModelOptions::set(util::MOLINEAUX_PAIRWISE_SAMPLE);
        }else{
            ETS_ASSERT( false );        // stop this test
        }
        if( repl_gamma ){
            ModelOptions::set(util::PARASITE_REPLICATION_GAMMA);
        }
        
        // Set parameters; all of these were estimated externally from OpenMalaria
        Parameters params(prepareParameters());
        
        // This sets up the model based on parameters and options
        WithinHost::MolineauxInfection::init(params);
    }
    
    static void MosqLifeCycle_init() {
        ModelOptions::reset();
        ModelOptions::set(util::VECTOR_LIFE_CYCLE_MODEL);
    }
    
    static double getPrescribedMg( const PkPd::LSTMModel& pkpd ){
        double r = 0.0;
        foreach( const PkPd::MedicateData& md, pkpd.medicateQueue ){
            r += md.qty;
        }
        return r;
    }
    
    static void medicate(PkPd::LSTMModel& pkpd, size_t typeIndex, double qty,
                         double time, double duration, double bodyMass){
        pkpd.medicateDrug(typeIndex, qty, time, duration, bodyMass);
    }
    
    static void clearMedicateQueue( PkPd::LSTMModel& pkpd ){
        pkpd.medicateQueue.clear();
    }
    
    static auto_ptr<Host::Human> createHuman(SimTime dateOfBirth){
        return auto_ptr<Host::Human>( new Host::Human(dateOfBirth) );
    }
    static void setHumanWH(Host::Human& human, WithinHost::WHInterface *wh){
        human.withinHostModel = wh;
    }
};

#endif
