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
#include "Clinical/CaseManagementCommon.h"
#include "Clinical/ESDecisionTree.h"	// needed for ESDecisionMap
#include "Clinical/parser.h"
#include "Pathogenesis/State.h"
#include "WithinHost/WithinHostModel.h"
#include "Monitoring/Survey.h"
#include "inputData.h"

#include <cassert>
#include <list>
#include <limits>
#include <boost/unordered_map.hpp>

class ESCaseManagementSuite;
class ESDecisionTreeSuite;

namespace OM { namespace Clinical {
    using WithinHost::WithinHostModel;
    using boost::unordered_map;

/// Auxilliary output from running case management
struct CMAuxOutput {
    enum Hospitalisation {
	NONE, IMMEDIATE, DELAYED
    };
    
    Hospitalisation hospitalisation;	///< was case hospitalised immediately, after a delay, or not at all?
    bool RDT_used;	///< was an RDT test used?
    bool microscopy_used;	///< was a microscopy test used?
};

/// Data used for a withinHostModel->medicate() call
struct MedicateData {
    MedicateData () :
	qty(numeric_limits< double >::signaling_NaN()),
	cost_qty(numeric_limits< double >::signaling_NaN()),
	time(numeric_limits< double >::signaling_NaN())
    {}
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
	abbrev & stream;
	qty & stream;
	cost_qty & stream;
	time & stream;
    }
    
    string abbrev;	/// Drug abbreviation
    double qty;		/// Quantity of drug prescribed (mg?)
    double cost_qty;	/// Effective quantity prescribed, with respect to costs
    double time;	/// Time to medicate at (days from start of timestep, may be >= 1 (not this timestep))
};

/** A final treatment schedule (after application of applicable modifiers). */
struct ESTreatmentSchedule {
    ESTreatmentSchedule (const scnXml::HSESTreatmentSchedule& sched);
    
    /// Multiply the quantity of each medication based on the value of this map.
    void multiplyQty (const parser::SymbolValueMap&, bool affectsCost, const string& errObj);
    /// Delay the time of each medication based on the value of this map, in hours.
    void delay (const parser::SymbolValueMap&, const string& errObj);
    /// Remove medications not in time range (in hours) described by this map.
    void selectTimeRange (const parser::SymbolRangeMap&, bool affectsCost, const string& errObj);
    
    /// Add medications into medicate queue
    inline void apply (list<MedicateData>& medicateQueue) const {
	for (vector<MedicateData>::const_iterator it = medications.begin(); it != medications.end(); ++it)
	    medicateQueue.push_back (*it);
    }
    
    private:
	/// Data for each medicate() call.
	vector<MedicateData> medications;
};

/** A set of all modified forms of a treatment schedule. Corresponds to a
 * treatment (base schedule + modifiers) given in the XML.
 * 
 * The reason we can't just have all schedules (all modified variants) in one
 * list, is because each set of modifications has its own mask.
 *****************************************************************************/
class ESTreatment {
    public:
	/** Construct from a base schedule and modifiers.
	 * 
	 * Also neede the decision value map. */
	ESTreatment (const ESDecisionValueMap& dvMap,
		     const scnXml::HSESTreatment& elt,
		     list<string>& required);
	~ESTreatment();
	
	/** Given an input ESDecisionValue, find a variant of the base treatment
	 * schedule.
	 *
	 * May return NULL (handled by ESDecisionMap::getTreatment()). */
	ESTreatmentSchedule* getSchedule (ESDecisionValue&) const;
	
    private:
	typedef unordered_map<ESDecisionValue,ESTreatmentSchedule*> Schedules;
	Schedules schedules;
	ESDecisionValue schedulesMask;
};

/** Decision trees representation, mapping inputs to a ESTreatmentSchedule pointer.
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
	 * @param cm XML element describing probabilistic decisions and treatments
	 * @param complicated Determines whether hard-coded decisions for the
	 * uncomplicated or complicated case are added. */
	void initialize (const ::scnXml::HSESCaseManagement& cm, bool complicated);
        
        /** Run decision tree to arrive at an outcome.
	 *
	 * @returns An outcome as a binary-or'd list of decision values. */
        ESDecisionValue determine (ESHostData& hostData) const;
	
	/** Given a decision-tree outcome, return a corresponding treatment
	 * schedule. Return-value should always point to an existing
	 * ESTreatmentSchedule object (shouldn't be deleted by the caller).
	 * 
	 * If the treatment decision is but not found, or found but a treatment
	 * schedule is not, an error is thrown. */
	ESTreatmentSchedule* getSchedule (ESDecisionValue outcome) const;
        
	/// Return one of CMAuxOutput::Hospitalisation's values.
	inline CMAuxOutput::Hospitalisation hospitalisation (ESDecisionValue outcome) const{
	    ESDecisionValue masked = outcome & hospitalisation_mask;
	    if( masked == hospitalisation_immediate )
		return CMAuxOutput::IMMEDIATE;
	    else if( masked == hospitalisation_delayed )
		return CMAuxOutput::DELAYED;
	    else
		return CMAuxOutput::NONE;
	}
	
	/// Return true if the input outcome indicates an RDT was used
	inline bool RDT_used (ESDecisionValue outcome) const{
	    return (outcome & test_mask) == test_RDT;
	}
	/// Return true if the input outcome indicates microscopy was used
	inline bool microscopy_used (ESDecisionValue outcome) const{
	    return (outcome & test_mask) == test_microscopy;
	}
	
    private:
	// All data here should be set by ESCaseManagement::init(); don't checkpoint.
	
	ESDecisionValueMap dvMap;
	
        // Currently we walk through all decisions, required or not
        vector<ESDecisionTree*> decisions;
	
	typedef unordered_map<ESDecisionValue,ESTreatment*> Treatments;
	Treatments treatments;
	// Used to mask ESDecisionValues before lookup in treatments:
	ESDecisionValue treatmentsMask;
	ESDecisionValue hospitalisation_mask,
				    hospitalisation_immediate,
				    hospitalisation_delayed;
	ESDecisionValue test_mask, test_RDT, test_microscopy;
	
	friend class ::ESCaseManagementSuite;	// unittests
	friend class ::ESDecisionTreeSuite;
};


/** Tracks clinical status (sickness), does case management for new events,
 * medicates treatment, determines patient recovery, death and sequelae.
 *****************************************************************************/
class ESCaseManagement : public CaseManagementCommon {
    public:
	static void init ();
	
	/** Load health system data from initial data or an intervention's data (both from XML).
	* (Re)loads all data affected by this healthSystem element. */
	static void setHealthSystem (const scnXml::HealthSystem& healthSystem);
	
	static void cleanup ();
	
	static void massDrugAdministration(list< OM::Clinical::MedicateData >& medicateQueue);
	
        /** Runs through case management decisions, selects treatments and
         * applies them to the passed medicateQueue.
         * 
         * Returns: some extra info (see CMAuxOutput definition). */
	static CMAuxOutput execute (list<MedicateData>& medicateQueue,
				    Pathogenesis::State pgState,
				    WithinHost::WithinHostModel& withinHostModel,
				    double ageYears,
				    Monitoring::AgeGroup ageGroup);
	
    private:
	
	//BEGIN Static parameters — set by init()
        static ESDecisionMap uncomplicated, complicated;
	
	/// MDA dosage info; null pointer if not provided
	static ESTreatmentSchedule* mdaDoses;
	//END
};

} }
#endif