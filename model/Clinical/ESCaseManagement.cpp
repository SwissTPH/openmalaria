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

#include "Clinical/ESCaseManagement.h"
#include "Clinical/CMDecisionTree.h"
#include "Clinical/EventScheduler.h"
#include "Monitoring/Survey.h"
#include "util/errors.h"
#include "util/ModelOptions.h"

#include <set>
#include <sstream>
#include <boost/format.hpp>
#include <boost/assign/std/vector.hpp> // for 'operator+=()'

namespace OM { namespace Clinical {
    using namespace OM::util;
    using namespace boost::assign;
    using namespace Monitoring;
    using boost::format;

// -----  ESTreatmentSchedule  -----

ESTreatmentSchedule::ESTreatmentSchedule (const scnXml::PKPDSchedule& sched) {
    const ::scnXml::PKPDSchedule::MedicateSequence& mSeq = sched.getMedicate();
    medications.resize (mSeq.size ());
    for (size_t j = 0; j < mSeq.size(); ++j) {
	medications[j].abbrev = mSeq[j].getDrug();
	medications[j].qty = mSeq[j].getMg();
	medications[j].cost_qty = medications[j].qty;	// effective quantity w.r.t. cost starts off equal to quantity
	medications[j].time = mSeq[j].getHour() / 24.0;	// convert from hours to days
	if( mSeq[j].getDuration().present() ){
	    if( !(mSeq[j].getDuration().get() > 0.0) ){
		throw util::xml_scenario_error( "duration of an IV dose must be some positive amount of time" );
	    }
	    medications[j].duration = mSeq[j].getDuration().get() / 24.0;
	}
    }
}


// -----  ESDecisionMap  -----

class ESDecisionMapProcessor {
    ESDecisionMap::TreeType treeType;
    ESDecisionValueMap& dvMap;
    
    // Collection of all decisions not yet added into decisions or culled.
    // Filled by constructor, emptied by process().
    typedef map<string,CMDecisionTree*> DecisionList;
    DecisionList pending;
    
    // Set of all tests required. Filled by process().
    set<CMDecisionTree*> required;
    
    // Add by key "decision", and return true if added
    bool addToPending( CMDecisionTree* d ){
	return pending.insert( make_pair(d->decision, d) ).second;
    }
    void addRequires (const string& x){
	DecisionList::const_iterator it = pending.find( x );
	if( it == pending.end() ) {
	    ostringstream msg;
	    msg << "ESCaseManagement: decision " << x << " required (for ";
            if( treeType == ESDecisionMap::MDA ){
                msg << "MDA";
            }else if( treeType == ESDecisionMap::Uncomplicated ){
                msg << "uncomplicated";
            }else if( treeType == ESDecisionMap::Complicated ){
                msg << "complicated";
            }else{
                assert (false);
            }
            msg << " tree)";
	    throw xml_scenario_error(msg.str());
	}
	if( required.insert( it->second ).second ){	// if newly inserted...
	    BOOST_FOREACH( const string& dep, it->second->depends )
		addRequires( dep );
	}	// else it was already considered
    }
    
    inline bool hasAllDependencies (const CMDecisionTree* decision, const set<string>& dependencies) {
	BOOST_FOREACH ( const string& n, decision->depends ) {
	    if (dependencies.count(n) == 0) {
		return false;
	    }
	}
	return true;
    }
    
public:
    /* Read from XML into a temporary list. */
    ESDecisionMapProcessor (ESDecisionValueMap& dvm,
                            const ::scnXml::HSESCaseManagement& xmlCM,
                            ESDecisionMap::TreeType tt) :
        treeType(tt), dvMap(dvm)
    {
	// Assemble a list of all tests we might need to add
	if ( treeType == ESDecisionMap::Uncomplicated ) {
	    addToPending (new CMDTCaseType (dvMap));	// this test only makes sense for UC case
	}
	addToPending (new ESDecisionParasiteTest (dvMap));	// optimised away if not used
	BOOST_FOREACH ( const ::scnXml::HSESDecision& xmlDc, xmlCM.getDecisions().getDecision() ) {
	    if( !addToPending (CMDecisionTree::create (dvMap, xmlDc)) )
		throw xml_scenario_error((format("Case management: decision %1% described twice") %xmlDc.getName()).str());
	}
    }
    
    /* Filter and output. */
    void process (ESDecisionMap::Decisions& decisions, list<string>& requiredOutputs) {
	BOOST_FOREACH( const string& req, requiredOutputs ){
	    addRequires( req );
	}
	
	decisions.reserve (pending.size());
	set<string> added;	// names of decisions added; used since it's faster to lookup decision names here than linearly in "decisions"
	for (DecisionList::iterator it = pending.begin(); ;) {
	    if (it == pending.end ()) {
		if (pending.empty ())	// good, we're done
		    break;
		// else: some elements had unresolved dependencies; this should already have been caught by addRequires though
		throw TRACED_EXCEPTION_DEFAULT("ESCaseManagement: didn't catch all dependencies");
	    }
	    //cout << "Considering " << (*it)->decision << " with " << (*it)->depends.size()<<" dependencies"<<endl;
	    if (required.count( it->second ) == 0) {
		//Note: optimising out unwanted decisions like this has two
		// possible sources of error, which should be checked:
		// 1. empty slots in the decisions list
		// 2. That required decisions are mistakenly optimised out.
#ifdef WITHOUT_BOINC	// only print in non-BOINC mode to reduce stderr.txt size on server; this isn't a crucial message
		if( it->first != "result" && it->first != "case" ) {	// these are built-in tests which may not be used
		    cerr << "Warning: ESCaseManagement: decision " << it->first << " is unused (for ";
                    if( treeType == ESDecisionMap::MDA ){
                        cerr << "MDA";
                    }else if( treeType == ESDecisionMap::Uncomplicated ){
                        cerr << "uncomplicated";
                    }else if( treeType == ESDecisionMap::Complicated ){
                        cerr << "complicated";
                    }else{
                        assert (false);
                    }
                    cerr << " tree)"<<endl;
		}
#endif
		delete it->second;
		pending.erase(it);
		it = pending.begin();	// after erase, we must restart from beginning
	    } else if (hasAllDependencies (it->second, added)) {
		decisions.push_back( it->second );
		added.insert (it->second->decision);
		pending.erase (it);
		it = pending.begin ();	// restart from beginning as above; also means we know if we reach the end there shouldn't be any elements left
	    } else
		++it;
	}
    }
};
inline ESDecisionValue treatmentGetValue (const ESDecisionValueMap::value_map_t& vmap, const string& value) {
    ESDecisionValueMap::value_map_t::const_iterator it = vmap.find (value);
    if (it == vmap.end())	// value is unknown
	throw xml_scenario_error((format("Treatment description given for treatment %1% which isn't an output of \"treatment\" decision") %value).str());
    return it->second;
}
void ESDecisionMap::initialize (const ::scnXml::HSESCaseManagement& xmlCM, TreeType treeType, bool reinitialise) {
    if( reinitialise ){
        // This function is also used to load a new health-system from intervention data; therefore clear old data:
        dvMap.clear();
        decisions.clear();
        treatments.clear();
    }else{
        if( !dvMap.empty() || !decisions.empty() || !treatments.empty() )
            // multiple MDA descriptions are probably the only thing that will cause this...
            throw util::unimplemented_exception( "multiple MDA descriptions using 1-day TS decision tree" );
    }
    
    // Construct processor & read from XML.
    // Fills dvMap (which must be done before evaluating treatments).
    ESDecisionMapProcessor processor( dvMap, xmlCM, treeType );
    
    list<string> required;	// list required decisions, to avoid optimising out
    
    if( treeType == Complicated ){
	// Register required "hospitalisation" output, before decisions are created
	// (by ESDecisionMapProcessor). Effects: no transmission in hospital, delay
	// simply adds one day to time before returning to treatment seeking.
	required.push_back( "hospitalisation" );
	vector<string> tfValues;
	tfValues += "none","immediate","delayed";
	hospitalisation_mask = dvMap.add_decision_values( "hospitalisation", tfValues );
	hospitalisation_immediate = dvMap.get( "hospitalisation", "immediate" );
	hospitalisation_delayed = dvMap.get( "hospitalisation", "delayed" );
    }
    
    // "test" decision: determines diagnostic usage and reported.
    required.push_back( "test" );
    test_mask = dvMap.getDecisionMask( "test" );
    test_RDT = dvMap.get( "test", "RDT" );
    test_microscopy = dvMap.get( "test", "microscopy" );
    
    // "result" decision: result of diagnostic (added automatically, but we need the keys)
    diagnostic_mask = dvMap.getDecisionMask( "result" );
    diagnostic_negative = dvMap.get( "result", "negative" );
    diagnostic_positive = dvMap.get( "result", "positive" );
    
    if( util::ModelOptions::option(util::NON_MALARIA_FEVERS) ){
        required.push_back( "AB_provider" );
        vector<string> outcomes;
        outcomes += "none","facility","informal";
        AB_provider_mask = dvMap.add_decision_values( "AB_provider", outcomes );
        AB_provider_facility = dvMap.get( "AB_provider", "facility" );
        AB_provider_informal = dvMap.get( "AB_provider", "informal" );
    }
    
    // Read treatments
    required.push_back( "treatment" );
    tuple< ESDecisionValue, ESDecisionValueMap::value_map_t > mask_vmap_pair = dvMap.getDecision("treatment");
    treatmentsMask = mask_vmap_pair.get<0>();
    
    const ESDecisionValueMap::value_map_t& treatmentCodes = mask_vmap_pair.get<1>();
    BOOST_FOREACH( const ::scnXml::HSESTreatment& treatment, xmlCM.getTreatments().getTreatment() ){
        ESDecisionValue key = treatmentGetValue( treatmentCodes, treatment.getName() );
	treatments.insert( key, new ESTreatment( dvMap, treatment, required ) );
    }
    
    // Filter and add decisions (must be done after reading treatments):
    processor.process( decisions, required );
}

ESDecisionMap::~ESDecisionMap () {
}

ESDecisionValue ESDecisionMap::determine (const OM::Clinical::CMHostData& hostData) const {
    ESDecisionValue outcomes;	// initialized to 0
    // At this point, decisions is ordered such that all dependencies should be
    // met if evaluated in order, so we just do that.
    for( Decisions::const_iterator it=decisions.begin(); it!=decisions.end(); ++it ){
    // Pass determine the outcomes of previous decisions, filtered to the decisions it depends on.
	// Get back another outcome and add it into the outcomes set.
	outcomes |= it->determine (outcomes, hostData);
    }
    return outcomes;
}
ESTreatmentSchedule& ESDecisionMap::getSchedule (ESDecisionValue outcome) {
    ESDecisionValue masked = outcome & treatmentsMask;
    Treatments::iterator it = treatments.find (masked);
    if (it != treatments.end ()) {
	ESTreatmentSchedule* ret = it->second->getSchedule( outcome );
	if( ret != NULL )
	    return *ret;
	throw TRACED_EXCEPTION_DEFAULT( "a required modifier decision's output is unexpected" );
    }
    
    ostringstream msg;
    msg<<"decision outcome "<<dvMap.format( masked )<<" not found in list of treatments";
    throw xml_scenario_error (msg.str());
}


// -----  ESCaseManagement  -----

ESDecisionMap ESCaseManagement::uncomplicated, ESCaseManagement::complicated;
ESDecisionMap ESCaseManagement::mda;
// ESDecisionMap decisionsUC, decisionsSev, decisionsMDA;

void ESCaseManagement::setHealthSystem( const scnXml::HealthSystem& healthSystem) {
    if( !healthSystem.getEventScheduler().present() )
	throw util::xml_scenario_error ("Expected EventScheduler section in healthSystem data (initial or intervention)");
    const scnXml::HSEventScheduler& esData = healthSystem.getEventScheduler().get();
    uncomplicated.initialize (esData.getUncomplicated (), ESDecisionMap::Uncomplicated, true);
    complicated.initialize (esData.getComplicated (), ESDecisionMap::Complicated, true);
    
    // Calling our parent class like this is messy. Changing this would require
    // moving change-of-health-system handling into ClinicalModel.
    ClinicalEventScheduler::setParameters( esData );
}

pair<ESDecisionValue, bool> executeTree(
        ESDecisionMap* map,
        const CMHostData& hostData,
        list<MedicateData>& medicateQueue
){
    map->execute (hostData);
    return make_pair( outcome, schedule.anyTreatments() );
}

void ESCaseManagement::initMDA (const scnXml::HSESCaseManagement& desc){
    mda.initialize( desc, ESDecisionMap::MDA, false );
}

void ESCaseManagement::massDrugAdministration(
        const CMHostData& hostData,
        list<MedicateData>& medicateQueue,
        const Host::Human& human,
        Monitoring::ReportMeasureI screeningReport,
        Monitoring::ReportMeasureI drugReport
){
    Survey::current().addInt( screeningReport, human, 1 );
    bool anyTreatment = executeTree( &mda, hostData, medicateQueue ).second;
    if( anyTreatment ){
        Survey::current().addInt( drugReport, human, 1 );
    }
}

CMAuxOutput ESCaseManagement::execute (
        const CMHostData& hostData,
        list<MedicateData>& medicateQueue
) {
    assert (hostData.pgState & Episode::SICK);
    // We always remove any queued medications.
    medicateQueue.clear();
    
    ESDecisionMap *map = (hostData.pgState & Episode::COMPLICATED) ? &complicated : &uncomplicated;
    ESDecisionValue outcome = executeTree( map, hostData, medicateQueue ).first;
    
    CMAuxOutput auxOut;
    auxOut.hospitalisation = CMAuxOutput::NONE;
    if( hostData.pgState & Episode::COMPLICATED )
	auxOut.hospitalisation = map->hospitalisation(outcome);
    auxOut.diagnostic = map->diagnostic(outcome);
    auxOut.AB_provider = map->AB_provider(outcome);
    return auxOut;
}

} }