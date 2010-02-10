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
#include "Clinical/ESDecision.h"
#include "Pathogenesis/State.h"
#include "WithinHost/WithinHostModel.h"
#include "Survey.h"
#include "inputData.h"

#include <cassert>
#include <list>
#include <boost/unordered_map.hpp>


namespace OM { namespace Clinical {

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
    CaseTreatment (const scnXml::CM_leaf::MedicateSequence mSeq) {
	medications.resize (mSeq.size ());
	for (size_t j = 0; j < mSeq.size(); ++j) {
	    medications[j].abbrev = mSeq[j].getName();
	    medications[j].qty = mSeq[j].getQty();
	    medications[j].time = mSeq[j].getTime() / 24.0;	// convert from hours to days
	}
    }
    
    /// Add medications into medicate queue
    inline void apply (list<MedicateData>& medicateQueue, cmid id) {
	// Extract treatment-seeking delay from id (branch of our case-management tree)
	int delay = (id & Decision::TSDELAY_MASK) >> Decision::TSDELAY_SHIFT;
	assert (delay <= Decision::TSDELAY_NUM_MAX);
	
	for (vector<MedicateData>::iterator it = medications.begin(); it != medications.end(); ++it) {
	    medicateQueue.push_back (*it);
	    medicateQueue.back().time += double(delay);
	}
    }
    
    /// Data for each medicate() call.
    vector<MedicateData> medications;
};

/// Pair of cmid and CaseTreatment&
// (would have used std::pair, but it can't store a reference
struct CaseTreatmentPair {
    CaseTreatmentPair (cmid id, CaseTreatment& ct) :
	first(id), second(ct)
    {}
    cmid first;
    CaseTreatment& second;
};

/** Tracks clinical status (sickness), does case management for new events,
 * medicates treatment, determines patient recovery, death and sequelae.
 */
class ESCaseManagement {
    public:
	static void init ();
	
	static void massDrugAdministration(list< OM::Clinical::MedicateData >& medicateQueue);
	
	static cmid execute (list<MedicateData>& medicateQueue, Pathogenesis::State pgState, WithinHost::WithinHostModel& withinHostModel, double ageYears, SurveyAgeGroup ageGroup);
	
    private:
	static CaseTreatmentPair traverse (cmid id);
	
	class CMNode {
	    public:
		virtual ~CMNode () {}
		virtual CaseTreatmentPair traverse (cmid id) =0;
	};
	class CMPBranchSet : public CMNode {
	    struct PBranch {
		cmid outcome;
		double cumP;
	    };
	    vector<PBranch> branches;	// must contain at least one entry; last must have cumP => 1.0
	    
	    public:
		CMPBranchSet (const scnXml::CM_pBranchSet::CM_pBranchSequence& branchSeq);
		
		virtual CaseTreatmentPair traverse (cmid id);
	};
	class CMLeaf : public CMNode {
	    CaseTreatment ct;
	    public:
		CMLeaf (CaseTreatment t) : ct(t) {}
		
		virtual CaseTreatmentPair traverse (cmid id);
	};
	
	//BEGIN Static parameters — set by init()
	typedef boost::unordered_map<cmid,CMNode*> TreeType;
	/// Tree probability-branches and leaf nodes.
	static TreeType cmTree;
	
	/// Mask applied to id before lookup in cmTree.
	static cmid cmMask;
	
	/// MDA dosage info; null pointer if not provided
	static CaseTreatment* mdaDoses;
	//END
};

} }
#endif