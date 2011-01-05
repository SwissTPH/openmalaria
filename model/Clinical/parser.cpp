/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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

#include "Clinical/parser.h"
#include "util/errors.h"

#include <sstream>

#include <boost/spirit/include/version.hpp>
#ifndef SPIRIT_VERSION
#error "No SPIRIT_VERSION macro (?)"
#elif SPIRIT_VERSION < 0x2010
#error "Spirit version 2.1 or later required!"
#endif

#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>

using namespace OM::util;
/** We use boost::spirit (2.1) for parsing here. A warning: debugging this
 * code when it won't compile is nearly impossible, so be very careful when
 * making changes!
 *************************************************************************/
namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

typedef std::pair<string,OM::Clinical::parser::DoubleRange> SymbolRangePair;
typedef std::vector<SymbolRangePair> SymbolRangeList;

// We need to tell fusion about our structs to make them first-class
// fusion citizens. This has to be in global scope.
BOOST_FUSION_ADAPT_STRUCT(
    OM::Clinical::parser::Branch,
    (OM::Clinical::parser::DecisionValue, dec_value)
    (OM::Clinical::parser::Outcome, outcome)
)
BOOST_FUSION_ADAPT_STRUCT(
    OM::Clinical::parser::BranchSet,
    (std::string, decision)
    (std::vector<OM::Clinical::parser::Branch>, branches)
)
BOOST_FUSION_ADAPT_STRUCT(
    OM::Clinical::parser::DoubleRange,
    (double, first)
    (double, second)
)
BOOST_FUSION_ADAPT_STRUCT(
    SymbolRangePair,
    (string, first)
    (OM::Clinical::parser::DoubleRange, second)
)

namespace OM { namespace Clinical {
    typedef string::const_iterator iter_t;
    namespace parser {
	template <typename Iterator>
	struct list_grammar : qi::grammar<Iterator, SymbolList(), ascii::space_type> {
	    list_grammar() : list_grammar::base_type(list, "value list") {
		using qi::alnum;
		using qi::lexeme;
		using phoenix::val;
		using namespace qi::labels;
		
		symbol %= lexeme[ +( alnum | '.' | '_' ) ];
		list %= symbol % ',';
		
		symbol.name( "symbol" );
		list.name( "list" );
	    }
	    
	    qi::rule<Iterator, string(), ascii::space_type> symbol;
	    qi::rule<Iterator, SymbolList(), ascii::space_type> list;
	};
	
	template <typename Iterator>
	struct DR_grammar : qi::grammar<Iterator, Outcome(), ascii::space_type>
	{
	    DR_grammar() : DR_grammar::base_type(tree, "decision tree")
	    {
		using qi::alnum;
		using qi::lexeme;
		using qi::lit;
		using qi::double_;
		using qi::eps;
		using phoenix::at_c;
		using phoenix::push_back;
		using phoenix::val;
		using phoenix::construct;
		using namespace qi::labels;
		
		symbol %= lexeme[ +( alnum | '.' | '_' ) ] ;	// symbol token
		value_symbol = symbol[ _val = _1 ];
		value_double = double_[ _val = _1 ];
		value_range = range[ _val = _1 ];
		range %= double_ > '-' > double_ ;
		tree %= eps >	// force a match: means empty input is an error, not no match
		    (branch_set | symbol )	// either a set of branches or a symbol
		;
		outcome %=			// as tree, but with braces or a colon
		    ('{' > tree > '}')
		    | (':' > symbol)
		;
		branch_p %= '(' > value_double > ')' > outcome ;
		branch_age %= '(' > value_range > ')' > outcome ;
		branch_decision %= '(' > value_symbol > ')' > outcome ;
		branch_set =
		    ( symbol >> &lit('(') )[	// match if we start with this
			at_c<0>(_val) = _1		// and set decision
		    ]
		    > (	// then we expect a branch
			( eps(at_c<0>(_val) == "p") > branch_p ) |
			( eps(at_c<0>(_val) == "age") > branch_age ) |
			branch_decision
		    )[
			push_back( at_c<1>(_val), _1 )
		    ]
		    > *(	// then, any number of times
			&symbol	// if we have a symbol
			> lit( at_c<0>(_val) )	// it must be decision
			> (	// followed by a branch
			    ( eps(at_c<0>(_val) == "p") > branch_p ) |
			    ( eps(at_c<0>(_val) == "age") > branch_age ) |
			    branch_decision
			)[
			    push_back( at_c<1>(_val), _1 )
			]
		    )
		;
		
		symbol.name( "symbol" );
		value_symbol.name( "value" );
		value_double.name( "value" );
		value_range.name( "value" );
		range.name( "range" );
		tree.name( "tree" );
		outcome.name( "outcome" );
		branch_p.name( "branch (p)" );
		branch_age.name( "branch (age)" );
		branch_decision.name( "branch (decision)" );
		branch_set.name( "branch_set" );
	    }
	    
	    qi::rule<Iterator, string(), ascii::space_type> symbol;
	    qi::rule<Iterator, DecisionValue(), ascii::space_type> value_symbol;
	    qi::rule<Iterator, DecisionValue(), ascii::space_type> value_double;
	    qi::rule<Iterator, DecisionValue(), ascii::space_type> value_range;
	    qi::rule<Iterator, DoubleRange(), ascii::space_type> range;
	    qi::rule<Iterator, Outcome(), ascii::space_type> tree;
	    qi::rule<Iterator, Outcome(), ascii::space_type> outcome;
	    qi::rule<Iterator, Branch(), ascii::space_type> branch_p;
	    qi::rule<Iterator, Branch(), ascii::space_type> branch_age;
	    qi::rule<Iterator, Branch(), ascii::space_type> branch_decision;
	    qi::rule<Iterator, BranchSet(), ascii::space_type> branch_set;
	};
	
	template <typename Iterator>
	struct SymbolValueMap_grammar : qi::grammar<Iterator, SymbolValueMap(), ascii::space_type> {
	    SymbolValueMap_grammar() : SymbolValueMap_grammar::base_type(map, "key-value map") {
		using qi::alnum;
		using qi::lexeme;
		using qi::double_;
		using qi::_val;
		using qi::_1;
		using qi::_2;
		using namespace phoenix;
		using namespace qi::labels;
		
		symbol %= lexeme[ +( alnum | '.' | '_' ) ];
		// conveniently enters stuff into a map, but doesn't seem to work when value-type is not fundamental:
		map = (symbol > '(' > double_ > ')')[_val[_1] = _2] % ',';
		
		symbol.name( "symbol" );
		map.name( "map" );
	    }
	    
	    qi::rule<Iterator, string(), ascii::space_type> symbol;
	    qi::rule<Iterator, SymbolValueMap(), ascii::space_type> map;
	};
	
	template <typename Iterator>
	struct SymbolRangeList_grammar : qi::grammar<Iterator, SymbolRangeList(), ascii::space_type> {
	    SymbolRangeList_grammar() : SymbolRangeList_grammar::base_type(map, "key-range map") {
		using qi::alnum;
		using qi::digit;
		using qi::lexeme;
		using qi::double_;
		using qi::_val;
		using qi::_1;
		using qi::_2;
		using phoenix::val;
		using namespace qi::labels;
		
		symbol %= lexeme[ +( alnum | '.' | '_' ) ];
		range %= double_ > '-' > double_;
		pair %= symbol > '(' > range > ')';
		map %= pair % ',';
		
		symbol.name( "symbol" );
		range.name( "range" );
		pair.name( "key-range" );
		map.name( "map" );
	    }
	    
	    qi::rule<Iterator, string(), ascii::space_type> symbol;
	    qi::rule<Iterator, DoubleRange(), ascii::space_type> range;
	    qi::rule<Iterator, SymbolRangePair(), ascii::space_type> pair;
	    qi::rule<Iterator, SymbolRangeList(), ascii::space_type> map;
	};
    }
    
    
    template<class RuleType, class ReturnType>
    ReturnType parseInternal( const string& s, const string& errMsg, const string& errObj ){
	RuleType rule;
	iter_t first = s.begin(); // we need a copy of the iterator, not a temporary
	ReturnType ret;
	
	try{
	    qi::phrase_parse(
		    first, s.end(),	// iterators
		    rule,			// rule
		    ascii::space,		// space skipper
		    ret			// output (type must match rule's attribute type)
	    );
	}catch( qi::expectation_failure<iter_t> e ){
	    ostringstream msg;
	    msg
		<< errMsg << errObj
		<< "; expecting: " << e.what_
		<< " here: \"" << string(e.first,s.end()) << '"'
	    ;
	     throw xml_scenario_error( msg.str() );
	}
	if (first != s.end ()) {
	    ostringstream msg;
	    msg
		<< errMsg << errObj
		<< "; remainder: " << string(first,s.end())
	    ;
	    throw xml_scenario_error (msg.str());
	}
	
	return ret;
    }
    
    parser::SymbolList parser::parseSymbolList (const string& s, const string& errObj) {
	return parseInternal<
	    parser::list_grammar<iter_t>,
	    SymbolList
	>(
	    s,
	    "failed to parse comma-separated fields for ",
	    errObj
	);
    }
    
    parser::Outcome parser::parseTree (const string& s, const string& errObj) {
	return parseInternal<
	    parser::DR_grammar<iter_t>,
	    parser::Outcome
	>(
	    s,
	    "failed to parse tree for ",
	    errObj
	);
    }
    
    parser::SymbolValueMap parser::parseSymbolValueMap (const string& s, const string& errObj) {
	return parseInternal<
	    parser::SymbolValueMap_grammar<iter_t>,
	    parser::SymbolValueMap
	>(
	    s,
	    "failed to parse comma-separated fields for ",
	    errObj
	);
    }
    
    parser::SymbolRangeMap parser::parseSymbolRangeMap (const string& s, const string& errObj) {
	SymbolRangeList list = parseInternal<
	    parser::SymbolRangeList_grammar<iter_t>,
	    SymbolRangeList
	>(
	    s,
	    "failed to parse comma-separated ranges for ",
	    errObj
	);
	
	// Now convert to a map (filling directly didn't work):
	SymbolRangeMap ret;
	for( SymbolRangeList::const_iterator it = list.begin(); it != list.end(); ++it )
	    ret[it->first] = it->second;
	
	return ret;
    }
} }
