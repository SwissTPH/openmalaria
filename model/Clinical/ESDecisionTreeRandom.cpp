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
#include <sstream>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

using namespace OM::util;

namespace OM { namespace Clinical {
    struct DR_processor {
	DR_processor (const string& n, ESDecisionValueMap& dvm, ESDecisionRandom& d) : name(n), dvMap(dvm), dR(d) {}
	void process (const parser::Outcome& outcome) {
	    processOutcome(outcome,
			   ESDecisionValue(),	/* 0: start with no outcomes */
			   1.0				/* probability of reaching here, given all required values */
	    );
	    checkProbabilities ();
	}
	
    private:
	void processBranches (const parser::Branches& branches,
			     ESDecisionValue dependValues,
			     double dependP
	){
	    assert (branches.size ());	// spirit parser shouldn't allow this to happen
	    string decision = branches[0].decision;
	    if (decision == "p") {
		
		double cum_p = 0.0;
		BOOST_FOREACH( const parser::Branch& branch, branches ) {
		    assert (branch.decision == decision);	//TODO: catch in parser
		    
		    double p = boost::lexical_cast<double>( branch.dec_value );
		    cum_p += p;
		    
		    processOutcome( branch.outcome, dependValues, dependP*p );
		}
		// Test cum_p is approx. 1.0 in case the input tree is wrong. In any case, we force probabilities to add to 1.0.
		if (cum_p < 0.999 || cum_p > 1.001)	//TODO: improve error message!
		    throw util::xml_scenario_error ("ESCaseManagement: probabilities don't add up to 1.0 for some set of probability branches");
		
	    } else {
		
		ESDecisionValueMap::value_map_t valMap = dvMap.getDecision( decision ).second;	// copy
		BOOST_FOREACH( const parser::Branch& branch, branches ) {
		    assert (branch.decision == decision);	//TODO: catch in parser
		    
		    ESDecisionValueMap::value_map_t::iterator valIt = valMap.find( branch.dec_value );
		    if( valIt == valMap.end() )
			throw xml_scenario_error((boost::format("%1%(%2%) encountered: %2% is not an outcome of %1%") %decision %branch.dec_value).str());
		    
		    processOutcome( branch.outcome, dependValues | valIt->second, dependP );
		    
		    valMap.erase( valIt );
		}
		if( !valMap.empty() ){	// error: not all options were included
		    //TODO: improve error message!
		    ostringstream msg;
		    msg << "Expected:";
		    for (ESDecisionValueMap::value_map_t::iterator it = valMap.begin(); it != valMap.end(); ++it)
			msg <<' '<<decision<<'('<<it->first<<')';
		    throw xml_scenario_error( msg.str() );
		}
		
	    }
	}
	void processOutcome (const parser::Outcome& outcome,
			     ESDecisionValue dependValues,
			     double dependP
	){
	    if ( const string* val_p = boost::get<string>( &outcome ) ) {
		ESDecisionValue val = dvMap.get( name, *val_p );
		size_t i = 0;	// get index i in dR.values of this outcome
		while (true) {
		    if (i >= dR.values.size())
			throw logic_error( ( boost::format("unable to find index for %1%(%2%) (code error)") %name %*val_p ).str() );
		    if( dR.values[i] == val )
			break;
		    ++i;
		}
		
		// find/make an entry for dependent decisions:
		//NOTE: valgrind complains about a memory leak _here_.. why?
		vector<double>& outcomes_cum_p = dR.map_cum_p[ dependValues ];
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
	    } else if ( const parser::Branches* brs_p = boost::get<parser::Branches>( &outcome ) ) {
		processBranches( *brs_p, dependValues, dependP );
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
	
	const string& name;	// name of tree this decision is for
	ESDecisionValueMap& dvMap;
	ESDecisionRandom& dR;
    };
    
    ESDecisionRandom::ESDecisionRandom (ESDecisionValueMap& dvMap, const ::scnXml::HSESDecision& xmlDc) {
	decision = xmlDc.getName();
	string decErrStr = "decision "+decision;
	
	depends = parser::parseSymbolList( xmlDc.getDepends(), decision+" depends attribute" );
	
	vector<string> valueList = parser::parseSymbolList( xmlDc.getValues(), decision+" values attribute" );
	// Calling a base-class function in the constructor, but it's not virtual so isn't an issue:
	setValues (dvMap, valueList);
	
	// Get content of this element:
	const ::xml_schema::String *content_p = dynamic_cast< const ::xml_schema::String * > (&xmlDc);
	if (content_p == NULL)
	    throw runtime_error ("ESDecision: bad upcast?!");
	
	DR_processor processor (decision, dvMap, *this);
	// 2-stage parse: first produces a parser::Outcome object, second the
	// processor object does the work of converting into map_cum_p:
	processor.process( parser::parseTree( *content_p, decision ) );
    }
    
} }