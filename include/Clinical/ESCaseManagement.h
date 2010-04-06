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
#include "Clinical/ESDecisionTree.h"	// needed for ESDecisionMap
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
    CaseTreatment (const ::scnXml::HSESTreatmentSchedule::MedicateSequence& mSeq) {
	medications.resize (mSeq.size ());
	for (size_t j = 0; j < mSeq.size(); ++j) {
	    medications[j].abbrev = mSeq[j].getDrug();
	    medications[j].qty = mSeq[j].getMg();
	    medications[j].time = mSeq[j].getHour() / 24.0;	// convert from hours to days
	}
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

/** Decision trees representation, mapping inputs to a CaseTreatment pointer.
 *
 * Used to represent a UC/UC2 or severe decision tree.
 *****************************************************************************/
class ESDecisionMap {
    public:
        /// Constructor. Element is created statically, so initialize later
        ESDecisionMap () {}
        ~ESDecisionMap();
	/** Read decision trees from an XML element.
	 *
	 * @param decisions XML element describing probabilistic decisions
	 * @param complicated Determines whether hard-coded decisions for the
	 * uncomplicated or complicated case are added. */
	void initialize (const ::scnXml::HSESCaseManagement& decisions, bool complicated);
        
        /** Run decision tree to arrive at an outcome.
	 *
	 * @returns An outcome as a binary-or'd list of decision values. */
        ESDecisionValue determine (ESHostData& hostData) const;
	
	/** Given a decision-tree outcome, return a corresponding case-treatment
	 * object. Return-value should always point to an existing CaseTreatment
	 * object (which shouldn't be deleted by the caller). */
	CaseTreatment* getTreatment (ESDecisionValue outcome) const;
        
    private:
	ESDecisionValueMap dvMap;
	
        // Currently we walk through all decisions, required or not
        vector<ESDecisionTree*> decisions;
	
	typedef unordered_map<ESDecisionValue,CaseTreatment*> treatments_t;
	treatments_t treatments;
	// Used to mask ESDecisionValues before lookup in treatments:
	ESDecisionValue treatmentMask;
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
        static ESDecisionMap uncomplicated, complicated;
	
	/// MDA dosage info; null pointer if not provided
	static CaseTreatment* mdaDoses;
	//END
};

} }
#endif