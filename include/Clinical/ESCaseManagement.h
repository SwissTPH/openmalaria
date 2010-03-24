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

#ifndef Hmod_ESCaseManagement
#define Hmod_ESCaseManagement

#include "Global.h"
#include "Pathogenesis/State.h"
#include "WithinHost/WithinHostModel.h"
#include "Survey.h"
#include "inputData.h"

#include <cassert>
#include <list>
#include <boost/unordered_map.hpp>


namespace OM { namespace Clinical {
    using WithinHost::WithinHostModel;
    using boost::unordered_map;

/// Data used for a withinHostModel->medicate() call
struct MedicateData {
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
	abbrev & stream;
	qty & stream;
	time & stream;
    }
    
    string abbrev;	/// Drug abbreviation
    double qty;		/// Quantity of drug prescribed (mg?)
    double time;	/// Time to medicate at (days from start of timestep, may be >= 1 (not this timestep))
};

/// Data type stored in decisions
struct CaseTreatment {
    CaseTreatment (/*const scnXml::CM_leaf::MedicateSequence mSeq*/) {
	//FIXME: get data
// 	medications.resize (mSeq.size ());
// 	for (size_t j = 0; j < mSeq.size(); ++j) {
// 	    medications[j].abbrev = mSeq[j].getName();
// 	    medications[j].qty = mSeq[j].getQty();
// 	    medications[j].time = mSeq[j].getTime() / 24.0;	// convert from hours to days
// 	}
    }
    
    /// Add medications into medicate queue
    inline void apply (list<MedicateData>& medicateQueue) {
	//TODO: apply treatment adjustments (missed doses, quality, age, delay (below))
	// Extract treatment-seeking delay from id (branch of our case-management tree)
// 	int delay = (id & Decision::TSDELAY_MASK) >> Decision::TSDELAY_SHIFT;
// 	assert (delay <= Decision::TSDELAY_NUM_MAX);
	
	for (vector<MedicateData>::iterator it = medications.begin(); it != medications.end(); ++it) {
	    medicateQueue.push_back (*it);
// 	    medicateQueue.back().time += double(delay);
	}
    }
    
    /// Data for each medicate() call.
    vector<MedicateData> medications;
};


/// Type used as identifier in decision trees.
typedef uint64_t ESDecisionID;
namespace Decision { enum DecisionEnums {
    /* Values here are written in hexadecimal: http://en.wikipedia.org/wiki/Hexadecimal
    * Many are designed to be "flags", so the value corresponds to a single bit:
    * http://en.wikipedia.org/wiki/Flag_byte
    * (note & | ^ are C++'s binary AND, OR and XOR operators).
    * Maximum value I'd recommend using as a flag: 0x4000_0000 */
    NONE                                = 0x0,
    
    /** MUST correspond to Pathogenesis::MORBIDITY_MASK.
    *
    * Pathogenesis constants from Constant.h within this mask may be used. */
    MORBIDITY_MASK              = 0x3F,
    
    AGE_UNDER5                  = 0x40,
    AGE_OVER5                   = 0x80,
    
    TEST_NONE                   = 0x0,
    TEST_MICROSCOPY     = 0x100,
    TEST_RDT                    = 0x200,
    TEST_MASK                   = 0x300,
    
    RESULT_POSITIVE     = 0x400,
    RESULT_NEGATIVE     = 0x800,
    RESULT_MASK         = 0xC00,
    
    TSDELAY_NUM_MAX     = 3,
    TSDELAY_SHIFT       = 24,
    TSDELAY_MASK        = 0x3000000,
}; }

struct ESDecisionName {
    ESDecisionName () : id(0) {}
    ESDecisionName (const char* name) {
        (*this) = name;
    }
    void operator= (const char* name);
    inline bool operator== (ESDecisionName& that) {
        return id == that.id;
    }
private:
    uint32_t id;
    static map<string,uint32_t> id_map;
    static uint32_t next_free;
};

struct ESHostData {
    ESHostData (double aY, WithinHostModel& wH) :
        ageYears(aY), withinHost(wH) {}
    double ageYears;
    WithinHostModel& withinHost;
};

/** Representation of one decision, random or deterministic (deterministic
 * decisions are hard-coded).
 *
 * Implementations (extending this base) are in the cpp file since they needn't
 * be shared.
 *****************************************************************************/
class ESDecisionTree {
    public:
        /// Run decision tree, with input filtered by mask.
        virtual ESDecisionID determine (ESDecisionID input, ESHostData& hostData) =0;
        
        vector<ESDecisionName> depends;      // other decisions this depends upon
        ESDecisionID mask;      // mask covering all outputs
        vector<ESDecisionID> values;    // ids associated with each possible output
};


/** Decision trees representation, mapping inputs to a CaseTreatment pointer.
 *
 * Used to represent a UC/UC2 or severe decision tree.
 *****************************************************************************/
class ESDecisionMap {
    public:
        /// Constructor. Element is created statically, so initialize later
        ESDecisionMap () {}
        ~ESDecisionMap();
	/** Initialization.
	 *
	 * @param complicated Determines whether hard-coded decisions for the
	 * uncomplicated or complicated case are added. */
	void initialize (const ::scnXml::HSESCMDecisions& decisions, bool complicated);
        
        /// Run decision tree. TODO: consider non-treatment outputs.
        CaseTreatment* determine (Pathogenesis::State pgState, ESHostData& hostData);
        
    private:
        // Currently we walk through all decisions, required or not
        vector<ESDecisionTree*> decisions;
        unordered_map<ESDecisionID,CaseTreatment*> treatments;
};


/** Tracks clinical status (sickness), does case management for new events,
 * medicates treatment, determines patient recovery, death and sequelae.
 *****************************************************************************/
class ESCaseManagement {
    public:
	static void init ();
	
	static void massDrugAdministration(list< OM::Clinical::MedicateData >& medicateQueue);
	
        /** Runs through case management decisions, selects treatments and
         * applies them to the passed medicateQueue. */
	static void execute (list<MedicateData>& medicateQueue, Pathogenesis::State pgState, WithinHost::WithinHostModel& withinHostModel, double ageYears, SurveyAgeGroup ageGroup);
	
    private:
	
	//BEGIN Static parameters — set by init()
        //FIXME: initialize
	static ESDecisionMap UC, severe;
	
	/// MDA dosage info; null pointer if not provided
	static CaseTreatment* mdaDoses;
	//END
};

} }
#endif