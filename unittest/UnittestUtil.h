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

#include "PkPd/PkPdModel.h"
// #include "PkPd/HoshenPkPdModel.h"
#include "PkPd/LSTMPkPdModel.h"
#include "PkPd/LSTMTreatments.h"
#include "WithinHost/Infection/Infection.h"
#include "WithinHost/WHFalciparum.h"
#include "Monitoring/Surveys.h"

#include "schema/pharmacology.h"
#include "schema/monitoring.h"

using namespace OM;
using namespace WithinHost;

using ::OM::util::ModelOptions;

namespace OM {
    namespace WithinHost {
        extern bool opt_common_whm;
    }
}

class UnittestUtil {
public:
    static void initTime(int daysPerStep){
        sim::init( daysPerStep, 90.0 /*max human age*/ );
        // we could just use zero, but we may spot more errors by using some weird number
        sim::sim_time = sim::fromYearsN(83.2591);
    }
    // Initialise surveys, to the minimum required not to crash
    static void initSurveys(){
        scnXml::OptionSet opts;
        scnXml::Surveys surveys( numeric_limits<double>::quiet_NaN() /* detection limit */);
        scnXml::AgeGroup ageGroups( 0.0 /* lower bound */ );
        scnXml::Monitoring mon( opts, surveys, ageGroups, "no monitoring" );
        Monitoring::Surveys.init( mon );
    }
    
    static void PkPdSuiteSetup (PkPd::PkPdModel::ActiveModel modelID) {
	TimeStep::init( 1, 90.0 );	// I think the drug model is always going to be used with an interval of 1 day.
	ModelOptions::reset();
        ModelOptions::set(util::INCLUDES_PK_PD);
	
	//Note: we fudge this call since it's not so easy to falsely initialize scenario element.
	//PkPdModel::init ();
	
	PkPd::PkPdModel::activeModel = modelID;
	if (modelID == PkPd::PkPdModel::LSTM_PKPD) {
            // Drugs
 	    scnXml::Allele allele ( 1.0 /* initial_frequency */, 3.45 /* max_killing_rate */, 0.6654 /* IC50 */, 2.5 /* slope */, "sensitive" /* name */ );
	    
	    scnXml::PD pd;
	    pd.getAllele().push_back (allele);
	    
	    scnXml::PK pk ( 0.006654 /* negligible_concentration */, 19.254 /* half_life */, 20.8 /* vol_dist */ );
	    
	    scnXml::PKPDDrug drug ( pd, pk, "MF" /* abbrev */ );
	    
            scnXml::Drugs drugs;
            drugs.getDrug().push_back (drug);
	    
	    PkPd::LSTMDrugType::init (drugs);
            
            // Treatments
            scnXml::PKPDSchedule sched1("sched1");
            sched1.getMedicate().push_back(
                scnXml::PKPDMedication("MF", 6 /*mg*/, 0 /*hour*/));
            
            scnXml::PKPDSchedule sched2("sched2");
            sched2.getMedicate().push_back(
                scnXml::PKPDMedication("MF", 2 /*mg*/, 0 /*hour*/));
            sched2.getMedicate().push_back(
                scnXml::PKPDMedication("MF", 5 /*mg*/, 12 /*hour*/));
            
            // a very basic dosage table, so that we can test it does what's expected
            scnXml::PKPDDosages dosage1("dosage1");
            dosage1.getAge().push_back(scnXml::PKPDDosageRange(0 /*age lb*/,1 /*mult*/));
            dosage1.getAge().push_back(scnXml::PKPDDosageRange(5 /*age lb*/,5 /*mult*/));
            
            scnXml::Treatments treatments;
            treatments.getSchedule().push_back(sched1);
            treatments.getSchedule().push_back(sched2);
            treatments.getDosages().push_back(dosage1);
            PkPd::LSTMTreatments::init(treatments);
	} /*else if (modelID == PkPd::PkPdModel::HOSHEN_PKPD) {
	    PkPd::ProteomeManager::init ();
	    PkPd::HoshenDrugType::init();
	} */else {
	    assert (false);
	}
    }
    
    // For when infection parameters shouldn't be used; enforce by setting to NaNs.
    static void Infection_init_NaN () {
	Infection::latentP = sim::zero();
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
    static void Infection_init_1day () {
        // don't set immunity vars: currently they're not used by unittest
        Infection_init_NaN();
        Infection::latentP = sim::fromDays(15);
    }
    
    static void DescriptiveInfection_init () {
	TimeStep::init( 5, 90.0 );
        ModelOptions::reset();
        ModelOptions::set(util::INCLUDES_PK_PD);
    }
    
    static void EmpiricalWHM_setup () {
	TimeStep::init( 1, 90.0 );
        ModelOptions::reset();
        ModelOptions::set(util::EMPIRICAL_WITHIN_HOST_MODEL);
        OM::WithinHost::opt_common_whm = true;
    }
    
    static void AgeGroupInterpolation_init() {
        TimeStep::init( 5, 90.0 );
    }
    
    static void MosqLifeCycle_init() {
        ModelOptions::reset();
        ModelOptions::set(util::VECTOR_LIFE_CYCLE_MODEL);
    }
    
    static double getPrescribedMg( const PkPd::LSTMPkPdModel& pkpd ){
        double r = 0.0;
        foreach( const PkPd::MedicateData& md, pkpd.medicateQueue ){
            r += md.qty;
        }
        return r;
    }
    
    static void medicate(PkPd::LSTMPkPdModel& pkpd, size_t typeIndex, double qty,
                         double time, double duration, double bodyMass){
        pkpd.medicateDrug(typeIndex, qty, time, duration, bodyMass);
    }
    
    static void clearMedicateQueue( PkPd::LSTMPkPdModel& pkpd ){
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
