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

// This is a small file to separate out ESDecisionRandom's parsers.

#include "Clinical/ESDecisionTree.h"
#include "Clinical/parser.h"
#include "util/errors.hpp"

#include <string>
#include <set>
#include <sstream>
#include <algorithm>	// find
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

using namespace OM::util;

namespace OM { namespace Clinical {
    /* Input processor to set mask and map_cum_p. */
    struct DR_processor {
	DR_processor (const ESDecisionValueMap& dvm, ESDecisionRandom& d) : dvMap(dvm), dR(d) {
	    BOOST_FOREACH( const string& dependency, dR.depends ){
		pair<ESDecisionValue,const ESDecisionValueMap::value_map_t&> decMap = dvMap.getDecision( dependency );
		dR.mask |= decMap.first;	// add this decision's inputs into the mask
		
		// decMap.second is converted to an ESDecisionValueSet:
		inputDependencies.push_back( make_pair( dependency, decMap.second ) );
	    }
	}
	
	void process (const parser::Outcome& outcome) {
	    processOutcome(outcome,
			   set<string>(),	/* empty set: no decisions yet processed */
			   ESDecisionValue(),	/* 0: start with no outcomes */
			   1.0				/* probability of reaching here, given all required values */
	    );
	    checkProbabilities ();
	}
	
    private:
	void processBranches (const parser::BranchSet& branchSet,
			      set<string> usedInputs,	// needs to take a copy
			      ESDecisionValue dependValues,
			     double dependP
	){
	    if (branchSet.decision == "p") {
		
		double cum_p = 0.0;
		BOOST_FOREACH( const parser::Branch& branch, branchSet.branches ) {
		    double p = boost::lexical_cast<double>( branch.dec_value );
		    cum_p += p;
		    
		    processOutcome( branch.outcome, usedInputs, dependValues, dependP*p );
		}
		// Test cum_p is approx. 1.0 in case the input tree is wrong. In any case, we force probabilities to add to 1.0.
		if (cum_p < 0.999 || cum_p > 1.001)
		    throw util::xml_scenario_error ( (
			boost::format("decision tree %1%: expected probability sum to be 1.0 but found %2%")
			%dR.decision
			%cum_p
		    ).str() );
		
	    } else {
		// This check is to make sure dependencies are listed (other
		// code depends on this). TODO: check decisions aren't listed
		// unnecessarily.
		if( find( dR.depends.begin(), dR.depends.end(), branchSet.decision ) == dR.depends.end() )
		    throw xml_scenario_error( (
			boost::format( "decision tree %1%: %2% not listed as a dependency" )
			%dR.decision
			%branchSet.decision
		    ).str() );
		usedInputs.insert( branchSet.decision );
		
		ESDecisionValueMap::value_map_t valMap = dvMap.getDecision( branchSet.decision ).second;	// copy
		BOOST_FOREACH( const parser::Branch& branch, branchSet.branches ) {
		    ESDecisionValueMap::value_map_t::iterator valIt = valMap.find( branch.dec_value );
		    if( valIt == valMap.end() )
			throw xml_scenario_error( (
			    boost::format("decision tree %3%: %1%(%2%) encountered: %2% is not an outcome of %1%")
			    %branchSet.decision
			    %branch.dec_value
			    %dR.decision
			).str());
		    
		    processOutcome( branch.outcome, usedInputs, dependValues | valIt->second, dependP );
		    
		    valMap.erase( valIt );
		}
		if( !valMap.empty() ){	// error: not all options were included
		    ostringstream msg;
		    msg << "decision tree "<<dR.decision<<": expected branches:";
		    for (ESDecisionValueMap::value_map_t::iterator it = valMap.begin(); it != valMap.end(); ++it)
			msg <<' '<<branchSet.decision<<'('<<it->first<<')';
		    throw xml_scenario_error( msg.str() );
		}
		
	    }
	}
	void processOutcome (const parser::Outcome& outcome,
			     const set<string>& usedInputs,	// doesn't edit
			     ESDecisionValue dependValues,
			     double dependP
	){
	    if ( const string* val_p = boost::get<string>( &outcome ) ) {
		ESDecisionValue val = dvMap.get( dR.decision, *val_p );
		size_t i = 0;	// get index i in dR.values of this outcome
		while (true) {
		    if (i >= dR.values.size())
			throw logic_error( (
			    boost::format("decision tree %1%: unable to find index for value %2% (code error)")
			    %dR.decision
			    %*val_p
			).str() );
		    if( dR.values[i] == val )
			break;
		    ++i;
		}
		
		// all input values which match the decisions we have made
		ESDecisionValueSet inputValues = dependValues;
		for( list< pair< string, ESDecisionValueSet > >::const_iterator it = inputDependencies.begin(); it != inputDependencies.end(); ++it ){
		    if( usedInputs.count( it->first ) == 0 )
			// This decion was not decided; look at all permutations by it
			inputValues |= it->second;
		}
		
		BOOST_FOREACH( ESDecisionValue inputValue, inputValues.values ){
		    // find/make an entry for dependent decisions:
		    //NOTE: valgrind complains about a memory leak _here_.. why?
		    vector<double>& outcomes_cum_p = dR.map_cum_p[ inputValue ];
		    /* print cum-prob-array (part 1):
		    if( outcomes_cum_p.empty() ) cout << "new ";
		    cout << "outcome cumulative probabilities for "<<dvMap.format( val )<<" (index "<<i<<"): ";
		    */
		    // make sure it's size is correct (will need resizing if just inserted):
		    outcomes_cum_p.resize( dR.values.size(), 0.0 );	// any new entries have p(0.0)
		    for (size_t j = i; j < outcomes_cum_p.size(); ++j)
			outcomes_cum_p[j] += dependP;
		    /* print cum-prob-array (part 2):
		    for (size_t j = 0; j < outcomes_cum_p.size(); ++j)
			cout <<" "<<outcomes_cum_p[j];
		    cout<<endl;
		    */
		}
	    } else if ( const parser::BranchSet* bs_p = boost::get<parser::BranchSet>( &outcome ) ) {
		processBranches( *bs_p, usedInputs, dependValues, dependP );
	    } else {
		assert (false);
	    }
	}
	
	void checkProbabilities () {
	    size_t l = dR.values.size();
	    size_t l1 = l - 1;
	    for( ESDecisionRandom::map_cum_p_t::iterator val_cum_p = dR.map_cum_p.begin(); val_cum_p != dR.map_cum_p.end(); ++val_cum_p ) {
		assert( val_cum_p->second.size() == l );
		// We force the last value to 1.0. It should be roughly that anyway due to
		// previous checks; it doesn't really matter if this gives allows a small error.
		val_cum_p->second[l1] = 1.0;
	    }
	}
	
	const ESDecisionValueMap& dvMap;
	ESDecisionRandom& dR;
	list< pair< string, ESDecisionValueSet > > inputDependencies;
    };
    
    ESDecisionRandom::ESDecisionRandom (ESDecisionValueMap& dvm, const ::scnXml::HSESDecision& xmlDc) /*: dvMap(dvm)*/ {
	decision = xmlDc.getName();
	string decErrStr = "decision "+decision;
	
	depends = parser::parseSymbolList( xmlDc.getDepends(), decision+" depends attribute" );
	
	vector<string> valueList = parser::parseSymbolList( xmlDc.getValues(), decision+" values attribute" );
	
	dvm.add_decision_values (decision, valueList);
	values.resize (valueList.size());
	for( size_t i = 0; i < valueList.size(); ++i ) {
	    values[i] = dvm.get (decision, valueList[i]);
	    //cout<<"added: "<<dvm.format(values[i])<<endl;
	}
	
	// Get content of this element:
	const ::xml_schema::String *content_p = dynamic_cast< const ::xml_schema::String * > (&xmlDc);
	if (content_p == NULL)
	    throw runtime_error ("ESDecision: bad upcast?!");
	
	DR_processor processor (dvm, *this);
	// 2-stage parse: first produces a parser::Outcome object, second the
	// processor object does the work of converting into map_cum_p:
	processor.process( parser::parseTree( *content_p, decision ) );
    }
    
} }