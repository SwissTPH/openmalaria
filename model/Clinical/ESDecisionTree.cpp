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

#include "Clinical/ESDecisionTree.h"
#include "util/random.h"
#include "util/errors.hpp"

#include <set>
#include <boost/format.hpp>
#include <boost/spirit/include/qi.hpp>

namespace OM { namespace Clinical {
    using namespace OM::util;
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;
    
    void ESDecisionTree::setValues (ESDecisionValueMap dvMap, const vector< string >& valueList) {
	mask = dvMap.add_decision_values (decision, valueList);
	values.resize (valueList.size());
	size_t i = 0;
	BOOST_FOREACH ( const string& value, valueList ) {
	    values[i] = dvMap.get (decision, value);
	}
    }
    ESDecisionRandom::ESDecisionRandom (ESDecisionValueMap dvMap, const ::scnXml::Decision& xmlDc) {
	decision = xmlDc.getName();
	
	// Set depends, values. Start by defining a rule matching a symbol:
	qi::rule<string::iterator, string(), ascii::space_type> symbol = qi::lexeme[+(qi::alnum | '.' | '_')];
	
	string s = xmlDc.getDepends();
	string::iterator first = s.begin(); // we need a copy of the iterator, not a temporary
	// Parse s into depends; note that the "attribute type" of the
	// expression must match the type of depends (a vector<string>):
	qi::phrase_parse(first, s.end(),
			    (symbol % ','),
			    ascii::space,
			    depends);
	
	vector<string> valueList;
	s = xmlDc.getValues();
	first = s.begin();
	// Same as above, for valueList:
	qi::phrase_parse(first, s.end(),
			    (symbol % ','),
			    ascii::space,
			    valueList);
	
	setValues (dvMap, valueList);
	
	//TODO: parse tree
	//BRANCH_SET := BRANCH+
	//BRANCH := DECISION '(' VALUE ')' ( ':' OUTCOME | '{' TREE '}' )
	//OUTCOME, DECISION, VALUE := SYMBOL
	qi::rule<string::iterator, ascii::space_type> tree;
	qi::rule<string::iterator, ascii::space_type> branch = symbol >> '(' > symbol > ')' > ( ':' > symbol | '{' > tree > '}' );
	tree = +branch | symbol;
	
	const ::xml_schema::String *content_p = dynamic_cast< const ::xml_schema::String * > (&xmlDc);
	if (content_p == NULL)
	    throw runtime_error ("ESDecision: bad upcast?!");
	s = *content_p;
	cout << "Got content: "<<s<<endl;
	first = s.begin();
	// For now, we ignore output and just test it wil pass the tree
	qi::phrase_parse(first, s.end(),
			    tree,
			    ascii::space);
	if (first != s.end ()) {
	    ostringstream msg;
	    msg << "ESDecision: failed to parse tree; remainder: " << string(first,s.end());
	    throw xml_scenario_error (msg.str());
	}
    }
    
    ESDecisionValue ESDecisionRandom::determine (const ESDecisionValue input, ESHostData& hostData) const {
	unordered_map<ESDecisionValue,vector<double> >::const_iterator it = set_cum_p.find (input);
	if (it == set_cum_p.end())
	    return ESDecisionValue();	// no decision
	double sample = random::uniform_01 ();
	size_t i = 0;
	while (it->second[i] < sample)
	    ++i;
	return values[i];
    }
    
    ESDecisionUC2Test::ESDecisionUC2Test (ESDecisionValueMap dvMap) {
	decision = "pathogenesisState";
	vector< string > valueList (2, "UC1");
	valueList[1] = "UC2";
	setValues (dvMap, valueList);
    }
    ESDecisionValue ESDecisionUC2Test::determine (const ESDecisionValue, ESHostData& hostData) const {
	assert (hostData.pgState & Pathogenesis::SICK && !(hostData.pgState & Pathogenesis::COMPLICATED));
	if (hostData.pgState & Pathogenesis::SECOND_CASE)	//TODO: check this actually gets set
	    return values[1];
	else
	    return values[0];
    }
    
    ESDecisionAge5Test::ESDecisionAge5Test (ESDecisionValueMap dvMap) {
	decision = "age5Test";
	vector< string > valueList (2, "under5");
	valueList[1] = "over5";
	setValues (dvMap, valueList);
    }
    ESDecisionValue ESDecisionAge5Test::determine (const ESDecisionValue input, ESHostData& hostData) const {
	if (hostData.ageYears >= 5.0)
	    return values[1];
	else
	    return values[0];
    }
    
    ESDecisionParasiteTest::ESDecisionParasiteTest (ESDecisionValueMap dvMap) {
	decision = "result";
	depends.resize (1, "test");	// Note: we check in add_decision_values() that this test has the outcomes we expect
	
	vector< string > valueList (2, "negative");
	valueList[1] = "positive";
	setValues (dvMap, valueList);
	
	test_none = dvMap.get("test","none");
	test_microscopy = dvMap.get("test", "microscopy");
	test_RDT = dvMap.get("test","RDT");
    }
    
    ESDecisionValue ESDecisionParasiteTest::determine (const ESDecisionValue input, ESHostData& hostData) const {
	if (input == test_microscopy) {
	    // if (hostData.withinHost.getDensity () > ...)
	    //FIXME: or values[0]
	    return values[1];
	} else if (input == test_RDT) {
	    //FIXME: or values[0]
	    return values[0];
	} else
	    return ESDecisionValue();	// 0, no decision
    }
    

std::size_t hash_value(ESDecisionValue const& b) {
    boost::hash<int> hasher;
    return hasher(b.id);
}
ESDecisionValue ESDecisionValueMap::add_decision_values (const string& decision, const std::vector< string > values) {
    if (decision == "test") {	// Simple way to check "test" decision has expected values
	set<string> expected;
	expected.insert ("none");
	expected.insert ("microscopy");
	expected.insert ("RDT");
	BOOST_FOREACH ( const string& value, values ) {
	    set<string>::iterator it = expected.find (value);
	    if (it == expected.end())
		throw xml_scenario_error ((boost::format("CaseManagement: \"test\" has unexpected outcome: %1%") % value).str());
	    else
		expected.erase (it);
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
    if (n_bits+next_bit>=(sizeof(next_bit)*8))	// (only valid on 8-bit-per-byte architectures)
	throw runtime_error ("ESDecisionValue design: insufficient bits");
    // Now we've got enough bits to represent all outcomes, starting at next_bit
    // Zero always means "missing value", so text starts at our first non-zero value:
    uint64_t next=(1<<next_bit), step;
    step=next;
    BOOST_FOREACH ( const string& value, values ) {
	string name = (boost::format("%1%(%2%)") %decision %value).str();
	bool success = id_map.insert (pair<string,uint64_t>(name, next)).second;
	if (!success) {
	    ostringstream msg;
	    msg <<"ESDecisionValue: value \""<<name<<"\" doesn't have a unique name";
	    throw runtime_error (msg.str());
	}
	next += step;
    }
    assert (next <= (1<<(n_bits+next_bit)));
    
    // Set mask so bits which are used by values are 1:
    ESDecisionValue mask;
    for (size_t i = 0; i < n_bits; ++i) {
	mask |= (1<<next_bit);
	++next_bit;
    }
    assert ((next-step & ~mask.id) == 0);		// last used value should be completely masked
    assert ((1<<next_bit & mask.id) == 0);		// this bit should be off the end of the mask
    return mask;
}
ESDecisionValue ESDecisionValueMap::get (const string& decision, const string& value) {
    string name = (boost::format("%1%(%2%)") %decision %value).str();
    map<string,uint64_t>::const_iterator it = id_map.find (name);
    if (it == id_map.end()) {
	throw runtime_error ("ESDecisionValue: name used before an entry was created for it");
    } else {
	ESDecisionValue ret;
	ret.id = it->second;
	return ret;
    }
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