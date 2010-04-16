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
#include <boost/assign/std/vector.hpp> // for 'operator+=()'

namespace OM { namespace Clinical {
    using namespace OM::util;
    using namespace boost::assign;
    
    ESDecisionValue ESDecisionRandom::determine (const ESDecisionValue input, const ESHostData& hostData) const {
	map_cum_p_t::const_iterator it = map_cum_p.find (input);
	if (it == map_cum_p.end()){
	    // All possible input combinations should be in map_cum_p
	    throw logic_error( "ESDecisionRandom: input combination not found in map (code error)" );
	}
	double sample = random::uniform_01 ();
	size_t i = 0;
	while (it->second[i] < sample)
	    ++i;
	return values[i];
    }
    
    ESDecisionUC2Test::ESDecisionUC2Test (ESDecisionValueMap& dvMap) {
	decision = "case";
	vector< string > valueList;
	valueList += "UC1", "UC2";
	dvMap.add_decision_values( decision, valueList );
	UC1 = dvMap.get( decision, "UC1" );
	UC2 = dvMap.get( decision, "UC2" );
    }
    ESDecisionValue ESDecisionUC2Test::determine (const ESDecisionValue, const ESHostData& hostData) const {
	assert (hostData.pgState & Pathogenesis::SICK && !(hostData.pgState & Pathogenesis::COMPLICATED));
	if (hostData.pgState & Pathogenesis::SECOND_CASE)
	    return UC2;
	else
	    return UC1;
    }
    
    ESDecisionAge5Test::ESDecisionAge5Test (ESDecisionValueMap& dvMap) {
	decision = "age5Test";
	vector< string > valueList;
	valueList += "under5", "over5";
	dvMap.add_decision_values( decision, valueList );
	under5 = dvMap.get( decision, "under5" );
	over5 = dvMap.get( decision, "over5" );
    }
    ESDecisionValue ESDecisionAge5Test::determine (const ESDecisionValue, const ESHostData& hostData) const {
	if (hostData.ageYears >= 5.0)
	    return over5;
	else
	    return under5;
    }
    
    ESDecisionParasiteTest::ESDecisionParasiteTest (ESDecisionValueMap& dvMap) {
	decision = "result";
	
	vector<string> testValues;
	testValues += "none", "microscopy", "RDT";
	// Add values, which (1) lets us create test_none, etc., (2) introduces a
	// check when the "test" decision is added later on output values, and
	// (3) allows us to set mask.
	mask = dvMap.add_decision_values ("test", testValues);
	test_none = dvMap.get("test","none");
	test_microscopy = dvMap.get("test", "microscopy");
	test_RDT = dvMap.get("test","RDT");
	
	depends.resize (1, "test");	// Note: we check in add_decision_values() that this test has the outcomes we expect
	
	vector< string > valueList;
	valueList += "none", "negative", "positive";
	dvMap.add_decision_values( decision, valueList );
	none = dvMap.get( decision, "none" );
	negative = dvMap.get( decision, "negative" );
	positive = dvMap.get( decision, "positive" );
    }
    
    ESDecisionValue ESDecisionParasiteTest::determine (const ESDecisionValue input, const ESHostData& hostData) const {
	if (input == test_none)	// no test
	    return none;
	else {
	    double dens = hostData.withinHost.getTotalDensity ();
	    double pPositive = 0.0;	// chance of a positive result
	    if (input == test_microscopy) {
		// Microscopy sensitivity/specificity data in Africa;
		// Source: expert opinion — Allan Schapira
		if (dens > 100.)
		    pPositive = .9;
		else if (dens > 0.)
		    pPositive = .75;
		else
		    pPositive = 1.0 - .75;	// specificity
	    } else {
		assert (input == test_RDT);
		// RDT sensitivity/specificity for Plasmodium falciparum in Africa
		// Source: Murray et al (Clinical Microbiological Reviews, Jan. 2008)
		if (dens > 500.) {
		    if (dens > 5000.)
			pPositive = .997;
		    else if (dens > 1000.)
			pPositive = .992;
		    else	// dens > 500.
			pPositive = .926;
		} else {
		    if (dens > 100.)
			pPositive = .892;
		    else if (dens > 0.)
			pPositive = .539;
		    else	// !(dens > 0.)
			pPositive = 1.0 - .942;	// specificity
		}
	    }
	    
	    if (random::uniform_01() < pPositive)
		return positive;
	    else
		return negative;
	}
	assert(false);	// should return before here
    }
    

std::size_t hash_value(ESDecisionValue const& b) {
  boost::hash<ESDecisionValue::id_type> hasher;
    return hasher(b.id);
}
ESDecisionValue ESDecisionValueMap::add_decision_values (const string& decision, const std::vector< string > values) {
    set<string> valueSet;
    BOOST_FOREACH( const string& v, values )
	if( !valueSet.insert( v ).second )
	    throw xml_scenario_error( (boost::format( "CaseManagement: decision %1%'s value %2% in value list twice!" ) %decision %v).str() );
    
    pair<id_map_type::iterator,bool> dec_pair = id_map.insert( make_pair (decision, make_pair ( ESDecisionValue(), value_map_t() ) ) );
    value_map_t& valMap = dec_pair.first->second.second;	// alias new map
    if (dec_pair.second) {	// new entry; fill it
	
	// got length l = values.size(); want minimal n such that: 2^n >= l
	// that is, n >= log_2 (l)
	// so n = ceil (log_2 (l))
	uint32_t n_bits = (uint32_t)std::ceil( log( double(values.size()) ) / log( 2.0 ) );
	if( n_bits + next_bit >= sizeof(next_bit) * 8 )	// (only valid on 8-bit-per-byte architectures)
	    throw runtime_error ("ESDecisionValue design: insufficient bits");
	
	// Now we've got enough bits to represent all outcomes, starting at next_bit
	// Zero always means "missing value", so text starts at our first non-zero value:
	id_type next=0, step=(1<<next_bit);
	BOOST_FOREACH ( const string& value, values ) {
// 	    cout<<"ESDecisionValue: "<<decision<<'('<<value<<"): "<<next<<endl;
	    valMap[value] = ESDecisionValue(next);
	    next += step;
	}
	next_bit += n_bits;
	assert (next <= (1u<<next_bit));
	
	// Set mask so bits which are used by values are 1:
	ESDecisionValue mask;
	for (value_map_t::const_iterator cur_val = valMap.begin(); cur_val != valMap.end(); ++cur_val)
	mask |= cur_val->second;
// 	cout<<"Mask for "<<decision<<": "<<mask<<endl;
	dec_pair.first->second.first = mask;
	
    } else {	// decision already exists; confirm values match
	
	set<string> new_values (values.begin(), values.end());
	for (value_map_t::const_iterator cur_val = valMap.begin(); cur_val != valMap.end(); ++cur_val) {
	    set<string>::iterator it = new_values.find (cur_val->first);
	    if (it == new_values.end())
		throw xml_scenario_error ((boost::format("CaseManagement: decision %1%'s values don't match; expected value: %2%") % decision % cur_val->first).str());
	    else
		new_values.erase (it);
	}
	if (!new_values.empty ()) {
	    ostringstream msg;
	    msg << "CaseManagement: decision "<<decision<<"'s values don't match; unexpected values:";
	    BOOST_FOREACH ( const string& v, new_values )
		msg << ' ' << v;
	    throw xml_scenario_error (msg.str());
	}
	
    }
    
    return dec_pair.first->second.first;
}
ESDecisionValue ESDecisionValueMap::get (const string& decision, const string& value) const {
    id_map_type::const_iterator it = id_map.find (decision);
    if (it == id_map.end())
	throw runtime_error ((boost::format("ESDecisionValueMap::get(): no decision %1%") %decision).str());
    
    value_map_t::const_iterator it2 = it->second.second.find (value);
    if (it2 == it->second.second.end())
	throw runtime_error ((boost::format("ESDecisionValueMap::get(): no value %1%(%2%)") %decision %value).str());
    
    //cout << "ESDecisionValueMap::get ("<<decision<<", "<<value<<"): "<<it2->second.id<<endl;
    return it2->second;
}
tuple< ESDecisionValue, const ESDecisionValueMap::value_map_t& > ESDecisionValueMap::getDecision (const string& decision) const {
    id_map_type::const_iterator it = id_map.find (decision);
    if (it == id_map.end ())
	throw invalid_argument ((boost::format ("ESDecisionValueMap: no decision %1%") %decision).str());
    return tuple<
      ESDecisionValue,
      const ESDecisionValueMap::value_map_t&
    >(it->second.first, it->second.second);
}

void ESDecisionValueMap::format( const ESDecisionValue v, ostream& stream ) const {
    bool second = false;	// prepend second, third, etc., with ", "
    for( id_map_type::const_iterator dec_it = id_map.begin(); dec_it != id_map.end(); ++dec_it ) {
	ESDecisionValue masked = v & dec_it->second.first;
	for( value_map_t::const_iterator it = dec_it->second.second.begin(); it != dec_it->second.second.end(); ++it ) {
	    if( masked == it->second ){
		if( second )
		    stream << ", ";
		stream << dec_it->first<<'('<<it->first<<')';
		second = true;
		goto foundValue;
	    }
	}
	assert( false );	// v matched mask but no value: this shouldn't happen!
	foundValue:;
    }
}

ESDecisionValueSet::ESDecisionValueSet (const ESDecisionValueMap::value_map_t valMap) {
    for( ESDecisionValueMap::value_map_t::const_iterator it = valMap.begin(); it != valMap.end(); ++it )
	values.push_back( it->second );
}

void ESDecisionValueSet::operator|= (const ESDecisionValueSet& that) {
    list<ESDecisionValue> old;
    old.swap( values );
    
    BOOST_FOREACH( ESDecisionValue v1, old ){
	BOOST_FOREACH( ESDecisionValue v2, that.values ){
	    values.push_back( v1 | v2 );
	}
    }
}

} }