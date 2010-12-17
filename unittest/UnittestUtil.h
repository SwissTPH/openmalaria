/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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
#include "WithinHost/Infection/Infection.h"
#include "WithinHost/WithinHostModel.h"

using namespace OM;
using namespace WithinHost;

class UnittestUtil {
public:
    static void PkPdSuiteSetup (PkPd::PkPdModel::ActiveModel modelID) {
	Global::interval = 1;	// I think the drug model is always going to be used with an interval of 1 day.
	util::ModelOptions::optSet = util::INCLUDES_PK_PD;
	
	//Note: we fudge this call since it's not so easy to falsely initialize scenario element.
	//PkPdModel::init ();
	
	PkPd::PkPdModel::activeModel = modelID;
        PkPd::PkPdModel::hetWeightMultStdDev = 0.0;
	if (modelID == PkPd::PkPdModel::LSTM_PKPD) {
            scnXml::AgeGroupValues agvElt;
            // We're not testing the interpolation, so a constant value is enough.
            // 60.0 would be the correct value (for age 21), but this is what
            // our old distribution gave us (avoids having to update results):
            agvElt.getGroup().push_back( scnXml::Group::Group( 55.4993, 0.0 ) );
            
            PkPd::PkPdModel::weight = util::AgeGroupInterpolation::makeObject( agvElt, "UnittestUtil_weight" );
            PkPd::PkPdModel::hetWeightMultStdDev = 0.0;
            // hetWeightMult must be large enough that birth weight is at least 0.5 kg:
            PkPd::PkPdModel::minHetWeightMult = 0.5 / (*PkPd::PkPdModel::weight)( 0.0 );
            
 	    scnXml::Allele allele ( 1.0 /* initial_frequency */, 3.45 /* max_killing_rate */, 0.6654 /* IC50 */, 2.5 /* slope */, "sensitive" /* name */ );
	    
	    scnXml::PD pd;
	    pd.getAllele().push_back (allele);
	    
	    scnXml::PK pk ( 0.006654 /* negligible_concentration */, 19.254 /* half_life */, 20.8 /* vol_dist */ );
	    
	    scnXml::Drug drug ( pd, pk, "MF" /* abbrev */ );
	    
	    scnXml::DrugDescription dd;
	    dd.getDrug().push_back (drug);
	    
	    PkPd::LSTMDrugType::init (dd);
	} else if (modelID == PkPd::PkPdModel::HOSHEN_PKPD) {
            assert( false );
// 	    PkPd::ProteomeManager::init ();
// 	    PkPd::HoshenDrugType::init();
	} else {
	    assert (false);
	}
    }
    static void PkPdSuiteTearDown () {
	if (PkPd::PkPdModel::activeModel == PkPd::PkPdModel::LSTM_PKPD) {
	    PkPd::LSTMDrugType::cleanup();
            util::AgeGroupInterpolation::freeObject( PkPd::PkPdModel::weight );
	} else if (PkPd::PkPdModel::activeModel == PkPd::PkPdModel::HOSHEN_PKPD) {
            assert( false );
// 	    PkPd::HoshenDrugType::cleanup();
// 	    PkPd::ProteomeManager::cleanup();
	}
    }
    
    // For when infection parameters shouldn't be used; enforce by setting to NaNs.
    static void Infection_init_NaN () {
	Infection::latentp = 0;
	Infection::cumulativeYstar = numeric_limits<float>::quiet_NaN();
	Infection::cumulativeHstar = numeric_limits<float>::quiet_NaN();
	Infection::alpha_m = numeric_limits<double>::quiet_NaN();
	Infection::decayM = numeric_limits<double>::quiet_NaN();
    }
    static void Infection_init () {
	// Note: these values were pulled from one source and shouldn't be taken as authoritative
	Infection::latentp = 3;
	Infection::cumulativeYstar = (float) 68564384.7102;
	Infection::cumulativeHstar = (float) 71.676733;
	Infection::alpha_m = 1.0 - exp(- 2.411434);
	Infection::decayM = 2.717773;
    }
    
    static void DescriptiveInfection_init () {
	Global::interval = 5;
	util::ModelOptions::optSet = util::INCLUDES_PK_PD;
    }
    
    static void EmpiricalWHM_setup () {
	Global::interval = 1;
	util::ModelOptions::optSet = util::EMPIRICAL_WITHIN_HOST_MODEL;
    }
    
    static void AgeGroupInterpolation_init() {
        Global::interval = 5;
        Global::intervalsPerYear = Global::DAYS_IN_YEAR/Global::interval;
        Global::yearsPerInterval = double(Global::interval) / double(Global::DAYS_IN_YEAR);
        Global::maxAgeIntervals = static_cast<int> (90.0 * Global::intervalsPerYear);
    }
    
    // only point of this function is that we give UnittestUtil "friend" status, not all unittest classes
    static void setTotalParasiteDensity (WithinHostModel& whm, double density) {
	whm.totalDensity = density;
    }
};

#endif