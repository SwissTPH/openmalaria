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
#include "inputData.h"
#include "util/random.h"
#include "util/errors.hpp"

#include <set>
#include <boost/format.hpp>
#include <boost/static_assert.hpp>
#include <boost/unordered_map.hpp>
#include <boost/spirit/include/qi.hpp>

namespace OM { namespace Clinical {
    using namespace OM::util;
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;

// -----  Types  -----

void ESDecisionTree::setValues (const vector< string >& valueList) {
    mask = ESDecisionValue::add_decision_values (decision, valueList);
    values.resize (valueList.size());
    size_t i = 0;
    BOOST_FOREACH ( const string& value, valueList ) {
	values[i].assign (decision, value);
    }
}
class ESDecisionRandom : public ESDecisionTree {
    public:
	ESDecisionRandom (const ::scnXml::Decision& xmlDc) {
	    decision = xmlDc.getName();
	    
	    // Set depends, values:
	    string s = xmlDc.getDepends();
	    string::iterator first = s.begin(); // we need a copy of the iterator, not a temporary
	    // Parse s into depends; note that the "attribute type" of the
	    // expression must match the type of depends (a vector<string>):
	    qi::phrase_parse(first, s.end(),
			     (qi::lexeme[+(qi::alnum | '.' | '_')] % ','),
			     ascii::space,
			     depends);
	    
	    vector<string> valueList;
	    s = xmlDc.getValues();
	    first = s.begin();
	    qi::phrase_parse(first, s.end(),
			     (qi::lexeme[+(qi::alnum | '.' | '_')] % ','),
			     ascii::space,
			     valueList);
	    
	    setValues (valueList);
	    
	    //TODO: parse tree
	}
	
	virtual ESDecisionValue determine (const ESDecisionValue input, ESHostData& hostData) const {
	    unordered_map<ESDecisionValue,vector<double> >::const_iterator it = set_cum_p.find (input);
	    if (it == set_cum_p.end())
		return ESDecisionValue();	// no decision
	    double sample = random::uniform_01 ();
	    size_t i = 0;
	    while (it->second[i] < sample)
		++i;
	    return values[i];
	}
	
    private:
	
	// a map of depended decision values to
	// cumulative probabilities associated with values
	//NOTE: be interesting to compare performance with std::map
	unordered_map<ESDecisionValue,vector<double> > set_cum_p;
};
class ESDecisionUC2Test : public ESDecisionTree {
    public:
	ESDecisionUC2Test () {
	    decision = "pathogenesisState";
	    vector< string > valueList (2, "UC1");
	    valueList[1] = "UC2";
	    setValues (valueList);
	}
	virtual ESDecisionValue determine (const ESDecisionValue, ESHostData& hostData) const {
	    assert (hostData.pgState & Pathogenesis::SICK && !(hostData.pgState & Pathogenesis::COMPLICATED));
	    if (hostData.pgState & Pathogenesis::SECOND_CASE)	//TODO: check this actually gets set
		return values[1];
	    else
		return values[0];
	}
};
class ESDecisionAge5Test : public ESDecisionTree {
    public:
	ESDecisionAge5Test () {
	    decision = "age5Test";
	    vector< string > valueList (2, "under5");
	    valueList[1] = "over5";
	    setValues (valueList);
	}
	virtual ESDecisionValue determine (const ESDecisionValue input, ESHostData& hostData) const {
	    if (hostData.ageYears >= 5.0)
		return values[1];
	    else
		return values[0];
	}
};
class ESDecisionParasiteTest : public ESDecisionTree {
    public:
	ESDecisionParasiteTest () {
	    decision = "result";
	    depends.resize (1, "test");	// Note: we check in add_decision_values() that this test has the outcomes we expect
	    
	    vector< string > valueList (2, "negative");
	    valueList[1] = "positive";
	    setValues (valueList);
	}
	
	virtual ESDecisionValue determine (const ESDecisionValue input, ESHostData& hostData) const {
	static const ESDecisionValue TEST_MICROSCOPY("test", "microscopy"),
		TEST_RDT("test","RDT"),
		TEST_NONE("test","none");
	    if (input == TEST_MICROSCOPY) {
		// if (hostData.withinHost.getDensity () > ...)
		//FIXME: or values[0]
		return values[1];
	    } else if (input == TEST_RDT) {
		//FIXME: or values[0]
		return values[0];
	    } else
		return ESDecisionValue();	// 0, no decision
	}
};


// ESCaseManagement::TreeType ESCaseManagement::cmTree;
// cmid ESCaseManagement::cmMask;

ESDecisionMap ESCaseManagement::UC, ESCaseManagement::severe;
CaseTreatment* ESCaseManagement::mdaDoses;

// -----  Static  -----

void ESCaseManagement::init () {
    // Assume EventScheduler data was checked present:
    const scnXml::HSEventScheduler& xmlESCM = InputData().getHealthSystem().getEventScheduler().get();
    UC.initialize (xmlESCM.getUC (), false);
    severe.initialize (xmlESCM.getSevere (), true);
    
    // MDA Intervention data
    const scnXml::Interventions::MDADescriptionOptional mdaDesc = InputData().getInterventions().getMDADescription();
    if (mdaDesc.present()) {
	mdaDoses = new CaseTreatment (/*FIXME mdaDesc.get().getMedicate()*/);
    } else
	mdaDoses = NULL;
}

inline bool hasAllDependencies (const ESDecisionTree* decision, const set<string>& dependencies) {
    BOOST_FOREACH ( const string& n, decision->depends ) {
	if (!dependencies.count(n))
	    return false;
    }
    return true;
}
void ESDecisionMap::initialize (const ::scnXml::HSESCMDecisions& xmlDcs, bool complicated) {
    // Assemble a list of all tests we need to add
    list<ESDecisionTree*> toAdd;
    toAdd.push_back (new ESDecisionAge5Test);
    if (!complicated) {
	toAdd.push_back (new ESDecisionUC2Test);
	toAdd.push_back (new ESDecisionParasiteTest);
    }
    BOOST_FOREACH ( const ::scnXml::Decision& xmlDc, xmlDcs.getDecision() ) {
	toAdd.push_back (new ESDecisionRandom (xmlDc));
    }
    
    decisions.resize (toAdd.size());
    size_t i = 0;
    set<string> added;	// names of decisions added; used since it's faster to lookup decision names here than linearly in "decisions"
    for (list<ESDecisionTree*>::iterator it = toAdd.begin(); ; ++it) {
	if (it == toAdd.end ()) {
	    if (toAdd.empty ())	// good, we're done
		break;
	    // else: some elements had unresolved dependencies
	    ostringstream msg;
	    msg << "ESCaseManagement: some decisions have unmeetable dependencies (for "<<(complicated?"":"un")<<"complicated tree):";
	    BOOST_FOREACH ( ESDecisionTree* decision, toAdd ) {
		msg << "decision: "<<decision->decision;
		msg << "\thas unmet dependency decisions:";
		BOOST_FOREACH ( string& dep, decision->depends ) {
		    if (!added.count(dep))
			msg << " " << dep;
		}
	    }
	    throw xml_scenario_error(msg.str());
	}
	if (hasAllDependencies (*it, added)) {
	    decisions[i++] = *it;
	    added.insert ((*it)->decision);
	    toAdd.erase (it);
	    it = toAdd.begin ();	// restart from beginning; affects order and means we know if we reach the end there shouldn't be any elements left
	}
    }
}
ESDecisionMap::~ESDecisionMap () {
    BOOST_FOREACH ( ESDecisionTree* d, decisions ) {
	delete d;
    }
}


void ESCaseManagement::massDrugAdministration(list<MedicateData>& medicateQueue) {
    if (mdaDoses == NULL)
	throw util::xml_scenario_error ("MDA intervention without description");
    else
	mdaDoses->apply(medicateQueue);
}

void ESDecisionName::operator= (const char* name) {
    map<string,uint32_t>::const_iterator it = id_map.find (name);
    if (it == id_map.end()) {
        id = next_free;
        next_free++;
        id_map[name] = id;
    } else {
        id = it->second;
    }
}
map<string,uint32_t> ESDecisionName::id_map;
uint32_t ESDecisionName::next_free = 1;

ESDecisionValue ESDecisionValue::add_decision_values (const string& decision, const std::vector< string > values) {
    if (decision == "test") {	// Simple way to check "test" decision has expected values
	set<string> expected;
	expected.insert ("none");
	expected.insert ("microscopy");
	expected.insert ("RDT");
	BOOST_FOREACH ( const string& value, values ) {
	    if (!expected.count (value))
		throw xml_scenario_error ((boost::format("CaseManagement: \"test\" has unexpected outcome: %1%") % value).str());
	}
	if (!expected.empty ()) {
	    ostringstream msg;
	    msg << "CaseManagement: expected \"test\" to have outcomes:";
	    BOOST_FOREACH ( const string& v, expected )
		msg << ' ' << v;
	    throw xml_scenario_error (msg.str());
	}
    }
    
    // got length l = values.size() + 1 (default, "no outcome"); want minimal n such that: 2^n >= l
    // that is, n >= log_2 (l)
    // so n = ceil (log_2 (l))
    uint32_t n_bits = std::ceil (log (values.size() + 1) / log(2.0));
    if (n_bits+next_bit>=(sizeof(next_bit)*8))
	throw runtime_error ("ESDecisionValue design: insufficient bits");
    uint64_t next=(1<<next_bit), step;
    step=next;
    BOOST_FOREACH ( const string& value, values ) {
	string name = (boost::format("%1%(%2%)") %decision %value).str();
	bool inserted = id_map.insert (pair<string,uint64_t>(name, next)).second;
	if (!inserted)
	    throw runtime_error ("ESDecisionValue: value doesn't have a unique name");
	next += step;
    }
    assert (next <= (1<<(n_bits+next_bit)));
    
    ESDecisionValue mask;
    for (size_t i = 0; i < n_bits; ++i) {
	mask.id |= (1<<next_bit);
	++next_bit;
    }
    assert ((next-step) & mask.id == 0);
    assert (next & mask.id != 0);
    return mask;
}
void ESDecisionValue::assign (const string& decision, const string& value) {
    string name = (boost::format("%1%(%2%)") %decision %value).str();
    map<string,uint64_t>::const_iterator it = id_map.find (name);
    if (it == id_map.end()) {
	throw runtime_error ("ESDecisionValue: name used before an entry was created for it");
    } else {
	id = it->second;
    }
}
std::size_t hash_value(ESDecisionValue const& b) {
    boost::hash<int> hasher;
    return hasher(b.id);
}
map<string,uint64_t> ESDecisionValue::id_map;
uint64_t ESDecisionValue::next_bit = 0;

CaseTreatment* ESDecisionMap::determine (OM::Clinical::ESHostData& hostData) {
    ESDecisionValue outcomes;	// initialized to 0
    // At this point, decisions is ordered such that all dependencies should be
    // met if evaluated in order, so we just do that.
    BOOST_FOREACH ( const ESDecisionTree* decision, decisions ) {
	// Pass determine the outcomes of previous decisions, filtered to the decisions it depends on.
	// Get back another outcome and add it into the outcomes set.
	outcomes |= decision->determine (outcomes & decision->mask, hostData);
    }
    // Find our outcome. FIXME: do we not want to mask this? It normally only depens on drug chosen, right?
    //TODO: suppositories as separate drug?
    unordered_map<ESDecisionValue,CaseTreatment*>::const_iterator treatment = treatments.find (outcomes);
    if (treatment == treatments.end ()) {
	//FIXME: what do we do? Complain or ignore?
    } else {
	return treatment->second;
    }
}


void ESCaseManagement::execute (list<MedicateData>& medicateQueue, Pathogenesis::State pgState, WithinHost::WithinHostModel& withinHostModel, double ageYears, SurveyAgeGroup ageGroup) {
    ESDecisionMap* map;
    assert (pgState & Pathogenesis::SICK);
    if (pgState & Pathogenesis::COMPLICATED)
        map = &severe;
    else
        map = &UC;
    
    ESHostData hostData (ageYears, withinHostModel, pgState);
    
    CaseTreatment* treatment = map->determine (hostData);
    
    // We always remove any queued medications.
    medicateQueue.clear();
    treatment->apply (medicateQueue);
}

// -----  CMNode derivatives  -----
/*
ESCaseManagement::CMPBranchSet::CMPBranchSet (const scnXml::CM_pBranchSet::CM_pBranchSequence& branchSeq) {
    double pAccumulation = 0.0;
    branches.resize (branchSeq.size());
    for (size_t i = 0; i < branchSeq.size(); ++i) {
	branches[i].outcome = branchSeq[i].getOutcome ();
	pAccumulation += branchSeq[i].getP ();
	branches[i].cumP = pAccumulation;
    }
    // Test cumP is approx. 1.0 (in case the XML is wrong).
    if (pAccumulation < 0.999 || pAccumulation > 1.001)
	throw util::xml_scenario_error ("EndPoint probabilities don't add up to 1.0 (CaseManagementTree)");
    // In any case, force it exactly 1.0 (because it could be slightly less,
    // meaning a random number x could have cumP<x<1.0, causing index errors.
    branches[branchSeq.size()-1].cumP = 1.0;
}*/

} }