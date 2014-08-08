/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2014 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2014 Liverpool School Of Tropical Medicine
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

#ifndef Hmod_ESCaseManagement
#define Hmod_ESCaseManagement

#include "Global.h"
#include "Clinical/CMDecisionTree.h"    // needed for ESDecisionMap
#include "WithinHost/WHInterface.h"
#include "Monitoring/Survey.h"
#include "schema/pharmacology.h"

#include <cassert>
#include <list>
#include <limits>
#include <boost/ptr_container/ptr_unordered_map.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

class ESCaseManagementSuite;
class ESDecisionTreeSuite;

namespace scnXml{
    
}
namespace OM { namespace Clinical {
    using WithinHost::WHInterface;
    using Monitoring::Survey;
    
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
    ESTreatmentSchedule (const scnXml::PKPDSchedule& sched);
    
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
        std::vector<MedicateData> medications;
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
         * uncomplicated or complicated case are added.
         * @param reinitialise If true, clear any previously initialised data
         * (e.g. for replacement health system); if false, throw an exception
         * if the map was previously initialised. */
        void initialize (const ::scnXml::DecisionTree& cm, TreeType treeType, bool reinitialise);
        
        /** Run decision tree to arrive at an outcome.
         *
         * @returns An outcome as a binary-or'd list of decision values. */
        void execute( CMHostData hostData ) const;
        
        typedef boost::ptr_vector<CMDecisionTree> Decisions;
    private:
        // All data here should be set by ESCaseManagement::init(); don't checkpoint.
        
//         ESDecisionValueMap dvMap;
        
        // Currently we walk through all decisions, required or not
        Decisions decisions;
        
//         typedef boost::ptr_unordered_map<ESDecisionValue,ESTreatment> Treatments;
//         Treatments treatments;
        
        friend class ::ESCaseManagementSuite;   // unittests
        friend class ::ESDecisionTreeSuite;
};


/** Tracks clinical status (sickness), does case management for new events,
 * medicates treatment, determines patient recovery, death and sequelae.
 *****************************************************************************/
class ESCaseManagement {
public:
    /** Load health system data from initial data or an intervention's data (both from XML).
    * (Re)loads all data affected by this healthSystem element. */
    static void setHealthSystem (const scnXml::HealthSystem& healthSystem);
    
    /** Set up MDA drug. Must be called if massDrugAdministration() is
        * ever used to deploy an MDA intervention. */
    static void initMDA (const scnXml::DecisionTree& desc);
    
    static void massDrugAdministration(
        const CMHostData& hostData,
        list<MedicateData>& medicateQueue,
        const Host::Human& human,
        Monitoring::ReportMeasureI screeningReport,
        Monitoring::ReportMeasureI drugReport
    );
    
    /** Runs through case management decisions, selects treatments and
        * applies them to the passed medicateQueue.
        * 
        * Returns: some extra info (see CMAuxOutput definition). */
    static CMAuxOutput execute (
        const CMHostData& hostData,
        list<MedicateData>& medicateQueue
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