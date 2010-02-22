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
#include "util/ModelOptions.hpp"

#include "PkPd/PkPdModel.h"
#include "PkPd/HoshenPkPdModel.h"
#include "PkPd/LSTMPkPdModel.h"
#include "WithinHost/Infection.h"

using namespace OM;
using WithinHost::Infection;

struct UnittestUtil {
    static void PkPdSuiteSetup (PkPd::PkPdModel::ActiveModel modelID) {
	Global::interval = 1;	// I think the drug model is always going to be used with an interval of 1 day.
	util::ModelOptions::optSet.set (util::INCLUDES_PK_PD);
	
	//Note: we fudge this call since it's not so easy to falsely initialize scenario element.
	//PkPdModel::init ();
	
	PkPd::PkPdModel::activeModel = modelID;
	if (modelID == PkPd::PkPdModel::LSTM_PKPD) {
	    scnXml::Allele allele ( 1.0 /* initial_frequency */, 3.45 /* max_killing_rate */, 0.6654 /* IC50 */, 2.5 /* slope */, "sensitive" /* name */ );
	    
	    scnXml::PD pd;
	    pd.getAllele().push_back (allele);
	    
	    scnXml::PK pk ( 0.006654 /* negligible_concentration */, 19.254 /* half_life */, 20.8 /* vol_dist */ );
	    
	    scnXml::Drug drug ( pd, pk, "MF" /* abbrev */ );
	    
	    scnXml::DrugDescription dd;
	    dd.getDrug().push_back (drug);
	    
	    PkPd::LSTMDrugType::init (dd);
	    PkPd::LSTMPkPdModel::init();
	} else if (modelID == PkPd::PkPdModel::HOSHEN_PKPD) {
	    PkPd::ProteomeManager::init ();
	    PkPd::HoshenDrugType::init();
	    PkPd::HoshenPkPdModel::init();
	} else {
	    assert (false);
	}
    }
    
    // For when infection parameters shouldn't be used; enforce by setting to NaNs.
    static void Infection_init_NaN () {
	Infection::cumulativeYstar = numeric_limits<double>::quiet_NaN();
	Infection::cumulativeHstar = numeric_limits<double>::quiet_NaN();
	Infection::alpha_m = numeric_limits<double>::quiet_NaN();
	Infection::decayM = numeric_limits<double>::quiet_NaN();
    }
};

#endif