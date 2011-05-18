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
#include "Clinical/ESDecisionTree.h"    // needed for ESDecisionMap
#include "Clinical/parser.h"
#include "Pathogenesis/State.h"
#include "WithinHost/WithinHostModel.h"
#include "Monitoring/Survey.h"

#include <cassert>
#include <list>
#include <limits>
#include <boost/ptr_container/ptr_unordered_map.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

class ESCaseManagementSuite;
class ESDecisionTreeSuite;

namespace scnXml{
    class HSESTreatmentSchedule;
    class HSESTreatment;
    class HSESCaseManagement;
}
namespace OM { namespace Clinical {
    using WithinHost::WithinHostModel;
    
/// Auxilliary output from running case management
struct CMAuxOutput {
    enum Hospitalisation {
        NONE=0, IMMEDIATE, DELAYED
    } hospitalisation;  ///< was case hospitalised immediately, after a delay, or not at all?
    enum Diagnostic {
        NO_TEST=0, POSITIVE, NEGATIVE
    } diagnostic;        ///< Was a malaria-parasite diagnostic used, and if so what was the outcome?
    enum AB_provider_T {
        NO_AB=0, FACILITY, INFORMAL
    } AB_provider;
};

/// Data used for a withinHostModel->medicate() call
struct MedicateData {
    MedicateData () :
        qty(numeric_limits< double >::signaling_NaN()),
        cost_qty(numeric_limits< double >::signaling_NaN()),
        time(numeric_limits< double >::signaling_NaN()),
        duration(numeric_limits< double >::quiet_NaN())
    {}
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        abbrev & stream;
        qty & stream;
        cost_qty & stream;
        time & stream;
        duration & stream;
    }
    
    string abbrev;      /// Drug abbreviation
    double qty;         /// Quantity of drug prescribed (mg?)
    double cost_qty;    /// Effective quantity prescribed, with respect to costs
    double time;        /// Time to medicate at (days from start of timestep, may be >= 1 (not this timestep))
    double duration;    /// Duration for IV purposes (use IV admin if a number, oral if is NaN)
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
    /// Does this contain a positive number of treatments?
    inline bool anyTreatments () const {
        return !medications.empty();
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
        ESTreatmentSchedule* getSchedule (ESDecisionValue&);
        
    private:
        typedef boost::ptr_unordered_map<ESDecisionValue,ESTreatmentSchedule> Schedules;
        Schedules schedules;
        ESDecisionValue schedulesMask;
};

/** Decision trees representation, mapping inputs to a ESTreatmentSchedule pointer.
 *
 * Used to represent a UC/UC2 or severe decision tree.
 *****************************************************************************/
class ESDecisionMap {
    public:
        enum TreeType {
            MDA,
            Uncomplicated,
            Complicated
        };
        
        /// Constructor. Element is created statically, so initialize later
        ESDecisionMap () {}
        ~ESDecisionMap();
        /** Read decision trees from an XML element.
         *
         * @param cm XML element describing probabilistic decisions and treatments
         * @param complicated Determines whether hard-coded decisions for the
         * uncomplicated or complicated case are added. */
        void initialize (const ::scnXml::HSESCaseManagement& cm, TreeType treeType);
        
        /** Run decision tree to arrive at an outcome.
         *
         * @returns An outcome as a binary-or'd list of decision values. */
        ESDecisionValue determine (const OM::Clinical::ESHostData& hostData) const;
        
        /** Given a decision-tree outcome, return a corresponding treatment
         * schedule.
         * 
         * If the treatment decision is not found, or found but a treatment
         * schedule is not, an error is thrown. */
        ESTreatmentSchedule& getSchedule (ESDecisionValue outcome);
        
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
        
        /// Return one of CMAuxOutput::Diagnostic's values.
        inline CMAuxOutput::Diagnostic diagnostic (ESDecisionValue outcome) const{
            ESDecisionValue masked = outcome & diagnostic_mask;
            if( masked == diagnostic_positive )
                return CMAuxOutput::POSITIVE;
            else if( masked == diagnostic_negative )
                return CMAuxOutput::NEGATIVE;
            else
                return CMAuxOutput::NO_TEST;
        }
        
        /// Return true if the input outcome indicates an RDT was used
        inline bool RDT_used (ESDecisionValue outcome) const{
            return (outcome & test_mask) == test_RDT;
        }
        /// Return true if the input outcome indicates microscopy was used
        inline bool microscopy_used (ESDecisionValue outcome) const{
            return (outcome & test_mask) == test_microscopy;
        }
        /** Return one of CMAuxOutput::AntibioticProvider's values. Will only
         * return a correct answer if NON_MALARIA_FEVERS option is enabled. */
        inline CMAuxOutput::AB_provider_T AB_provider (ESDecisionValue outcome) const{
            ESDecisionValue masked = outcome & AB_provider_mask;
            if( masked == AB_provider_facility )
                return CMAuxOutput::FACILITY;
            else if( masked == AB_provider_informal )
                return CMAuxOutput::INFORMAL;
            else
                return CMAuxOutput::NO_AB;
        }
        
        typedef boost::ptr_vector<ESDecisionTree> Decisions;
    private:
        // All data here should be set by ESCaseManagement::init(); don't checkpoint.
        
        ESDecisionValueMap dvMap;
        
        // Currently we walk through all decisions, required or not
        Decisions decisions;
        
        typedef boost::ptr_unordered_map<ESDecisionValue,ESTreatment> Treatments;
        Treatments treatments;
        // Used to mask ESDecisionValues before lookup in treatments:
        ESDecisionValue treatmentsMask;
        ESDecisionValue hospitalisation_mask,
                                    hospitalisation_immediate,
                                    hospitalisation_delayed;
        ESDecisionValue test_mask, test_RDT, test_microscopy;
        ESDecisionValue diagnostic_mask,
                                    diagnostic_positive,
                                    diagnostic_negative;
        ESDecisionValue AB_provider_mask,
                                    AB_provider_facility,
                                    AB_provider_informal;
        
        friend class ::ESCaseManagementSuite;   // unittests
        friend class ::ESDecisionTreeSuite;
};


/** Tracks clinical status (sickness), does case management for new events,
 * medicates treatment, determines patient recovery, death and sequelae.
 *****************************************************************************/
class ESCaseManagement : public CaseManagementCommon {
    public:
        /** Load health system data from initial data or an intervention's data (both from XML).
        * (Re)loads all data affected by this healthSystem element. */
        static void setHealthSystem (const scnXml::HealthSystem& healthSystem);
        
        /** Set up MDA drug. Must be called if massDrugAdministration() is
         * ever used to deploy an MDA intervention. */
        static void initMDA (const scnXml::MDA::DescriptionType& desc);
        
        static void massDrugAdministration(
            const ESHostData& hostData,
            list<MedicateData>& medicateQueue,
            bool inCohort,
            Monitoring::AgeGroup ageGroup
        );
        
        /** Runs through case management decisions, selects treatments and
         * applies them to the passed medicateQueue.
         * 
         * Returns: some extra info (see CMAuxOutput definition). */
        static CMAuxOutput execute (
            const ESHostData& hostData,
            list<MedicateData>& medicateQueue,
            bool inCohort
        );
        
    private:
        
        //BEGIN Static parameters — set by init()
        static ESDecisionMap uncomplicated, complicated;
        
        /// MDA description
        static ESDecisionMap mda;
        //END
};

} }
#endif