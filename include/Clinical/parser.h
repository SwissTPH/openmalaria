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

#ifndef Hmod_Clinical_parser
#define Hmod_Clinical_parser

#include <string>
#include <vector>
#include <map>
#include <boost/config/warning_disable.hpp>
#include <boost/variant.hpp>

using std::string;
using std::pair;

namespace OM { namespace Clinical {
    namespace parser {
	typedef std::vector<string> SymbolList;
	
	class Branch;
	typedef std::vector<Branch> Branches;
	typedef boost::variant< boost::recursive_wrapper< Branches >, string > Outcome;
	struct Branch {
	    string decision;
	    string dec_value;
	    Outcome outcome;
	};
	
	typedef std::map<string,double> SymbolValueMap;
	
	typedef std::pair<double,double> DoubleRange;
	typedef std::map<string,DoubleRange> SymbolRangeMap;
	
	
	/** @brief Parser functions.
	 *
	 * Each takes a string s which it parses, and a string errObj which it
	 * uses to explain where a parse error occurred. */
	//@{
	/** Parse s as a comma-separated list of symbols and return a
	 * SymbolList (a symbol must match [a-zA-Z0-9\._]+). */
	SymbolList parseSymbolList (const string& s, const string& errObj);
	
	/** Parse s as a decision tree. */
	Outcome parseTree (const string& s, const string& errObj);
	
	/** Parse s as a comma-separated list of symbol-value pairs of the form
	 * SYMBOL(VALUE), where symbol is as above and value is a number. */
	SymbolValueMap parseSymbolValueMap (const string& s, const string& errObj);
	
	/** Similar to parseSymbolValueMap, but value is a double pair of the
	 * form A-B. */
	SymbolRangeMap parseSymbolRangeMap (const string& s, const string& errObj);
	//@}
    }
} }
#endif