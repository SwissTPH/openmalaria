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

#include "Clinical/ESCaseManagement.h"
#include "Clinical/ESDecisionTree.h"
#include "inputData.h"
#include "util/errors.hpp"

#include <set>
#include <boost/format.hpp>

namespace OM { namespace Clinical {
    using namespace OM::util;
    using boost::format;


// -----  ESDecisionMap  -----

inline bool hasAllDependencies (const ESDecisionTree* decision, const set<string>& dependencies) {
    BOOST_FOREACH ( const string& n, decision->depends ) {
	if (dependencies.count(n) == 0) {
	    return false;
	}
    }
    return true;
}
inline ESDecisionValue treatmentGetValue (const ESDecisionValueMap::value_map_t& vmap, const string& value) {
    ESDecisionValueMap::value_map_t::const_iterator it = vmap.find (value);
    if (it == vmap.end())
	// value is "void" or something unknown; neither is acceptable
	throw xml_scenario_error((format("Treatment for drug %1% which isn't an output of \"drug\" decision") %value).str());
    return it->second;
}
void ESDecisionMap::initialize (const ::scnXml::HSESCaseManagement& xmlCM, bool complicated) {
    // Assemble a list of all tests we need to add
    list<ESDecisionTree*> toAdd;
    toAdd.push_back (new ESDecisionAge5Test (dvMap));
    if (!complicated) {
	toAdd.push_back (new ESDecisionUC2Test (dvMap));
	toAdd.push_back (new ESDecisionParasiteTest (dvMap));
    }
    BOOST_FOREACH ( const ::scnXml::Decision& xmlDc, xmlCM.getDecision() ) {
	toAdd.push_back (new ESDecisionRandom (dvMap, xmlDc));
    }
    
    decisions.resize (toAdd.size());
    size_t i = 0;
    set<string> added;	// names of decisions added; used since it's faster to lookup decision names here than linearly in "decisions"
    for (list<ESDecisionTree*>::iterator it = toAdd.begin(); ;) {
	if (it == toAdd.end ()) {
	    if (toAdd.empty ())	// good, we're done
		break;
	    // else: some elements had unresolved dependencies
	    cerr << "ESCaseManagement: some decisions have unmeetable dependencies (for "<<(complicated?"":"un")<<"complicated tree):";
	    BOOST_FOREACH ( ESDecisionTree* decision, toAdd ) {
		cerr << "\ndecision: "<<decision->decision;
		cerr << "\thas unmet dependency decisions:";
		BOOST_FOREACH ( string& dep, decision->depends ) {
		    if (added.count(dep))
			cerr << " (" << dep << ')';
		    else
			cerr << " " << dep;
		}
	    }
	    cerr << endl;
	    throw xml_scenario_error("ESCaseManagement: some decisions have unmeetable dependencies (see above)");
	}
	//cout << "Considering " << (*it)->decision << " with " << (*it)->depends.size()<<" dependencies"<<endl;
	if (hasAllDependencies (*it, added)) {
	    decisions[i++] = *it;
	    added.insert ((*it)->decision);
	    toAdd.erase (it);
	    it = toAdd.begin ();	// restart from beginning; affects order and means we know if we reach the end there shouldn't be any elements left
	} else
	    ++it;
    }
    
    
    // Read treatments:
    pair< ESDecisionValue, ESDecisionValueMap::value_map_t > mask_vmap_pair = dvMap.getDecision("drug");
    const ESDecisionValueMap::value_map_t& drugCodes = mask_vmap_pair.second;
    treatments_t unmodified_treatments;
    BOOST_FOREACH( const ::scnXml::HSESTreatment& treatment, xmlCM.getTreatments().getTreatment() ){
	unmodified_treatments.insert( make_pair( treatmentGetValue( drugCodes, treatment.getDrug() ), new CaseTreatment( treatment.getMedicate() ) ) );
    }
    //TODO: introduce and apply modifiers
    // for now, just copy
    treatments = unmodified_treatments;
    // Include "void" input with an empty CaseTreatment:
    treatments[ESDecisionValue()] = new CaseTreatment (::scnXml::HSESTreatment::MedicateSequence());
    
    treatmentMask = mask_vmap_pair.first;
}
ESDecisionMap::~ESDecisionMap () {
    BOOST_FOREACH ( ESDecisionTree* d, decisions ) {
	delete d;
    }
    for( treatments_t::iterator it = treatments.begin(); it != treatments.end(); ++it )
	delete it->second;
}

ESDecisionValue ESDecisionMap::determine (OM::Clinical::ESHostData& hostData) const {
    ESDecisionValue outcomes;	// initialized to 0
    // At this point, decisions is ordered such that all dependencies should be
    // met if evaluated in order, so we just do that.
    BOOST_FOREACH ( const ESDecisionTree* decision, decisions ) {
	// Pass determine the outcomes of previous decisions, filtered to the decisions it depends on.
	// Get back another outcome and add it into the outcomes set.
	outcomes |= decision->determine (outcomes & decision->mask, hostData);
    }
    return outcomes;
}
CaseTreatment* ESDecisionMap::getTreatment (ESDecisionValue outcome) const {
    // Find our outcome.
    //TODO: suppositories as separate drug?
    ESDecisionValue masked = outcome & treatmentMask;
    unordered_map<ESDecisionValue,CaseTreatment*>::const_iterator treatment = treatments.find (masked);
    if (treatment == treatments.end ()) {
	ostringstream msg;
	msg<<"decision outcome "<<dvMap.format( masked )<<" not found in list of treatments";
	throw xml_scenario_error (msg.str());
    } else {
	return treatment->second;
    }
}


// -----  ESCaseManagement  -----

ESDecisionMap ESCaseManagement::uncomplicated, ESCaseManagement::complicated;
CaseTreatment* ESCaseManagement::mdaDoses;

void ESCaseManagement::init () {
    // Assume EventScheduler data was checked present:
    const scnXml::HSEventScheduler& xmlESCM = InputData().getHealthSystem().getEventScheduler().get();
    uncomplicated.initialize (xmlESCM.getUncomplicated (), false);
    complicated.initialize (xmlESCM.getComplicated (), true);
    
    // MDA Intervention data
    const scnXml::Interventions::MDADescriptionOptional mdaDesc = InputData().getInterventions().getMDADescription();
    if (mdaDesc.present()) {
	mdaDoses = new CaseTreatment ( mdaDesc.get().getMedicate() );
    } else
	mdaDoses = NULL;
}

void ESCaseManagement::massDrugAdministration(list<MedicateData>& medicateQueue) {
    if (mdaDoses == NULL)
	throw util::xml_scenario_error ("MDA intervention without description");
    else
	mdaDoses->apply(medicateQueue);
}

void ESCaseManagement::execute (list<MedicateData>& medicateQueue, Pathogenesis::State pgState, WithinHost::WithinHostModel& withinHostModel, double ageYears, SurveyAgeGroup ageGroup) {
    ESDecisionMap* map;
    assert (pgState & Pathogenesis::SICK);
    if (pgState & Pathogenesis::COMPLICATED)
        map = &complicated;
    else
        map = &uncomplicated;
    
    ESHostData hostData (ageYears, withinHostModel, pgState);
    
    ESDecisionValue outcome = map->determine (hostData);
    
    CaseTreatment* treatment = map->getTreatment(outcome);
    
    // We always remove any queued medications.
    medicateQueue.clear();
    treatment->apply (medicateQueue);
}

} }