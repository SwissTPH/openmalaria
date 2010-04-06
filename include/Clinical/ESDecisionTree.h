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

#ifndef Hmod_ESDecisionTree
#define Hmod_ESDecisionTree

#include "Global.h"
#include "Pathogenesis/State.h"
#include "WithinHost/WithinHostModel.h"
#include "inputData.h"

#include <map>
#include <boost/unordered_map.hpp>

class ESDecisionTreeSuite;

namespace OM { namespace Clinical {
    using WithinHost::WithinHostModel;
    using boost::unordered_map;

/** A compressed representation of all decision outcomes.
 * Pass by value (it is 64-bits in size). */
struct ESDecisionValue {
    ESDecisionValue () : id(0) {}
    
    inline bool operator== (const ESDecisionValue that) const {
	return id == that.id;
    }
    inline ESDecisionValue operator& (const ESDecisionValue that) const {
	return ESDecisionValue(id & that.id);
    }
    inline ESDecisionValue operator| (const ESDecisionValue that) {
	return ESDecisionValue( id | that.id );
    }
    inline void operator|= (const ESDecisionValue that) {
	id |= that.id;
    }
    private:
	typedef uint64_t id_type;
	
	// private constructor, only for use by internal operations
	ESDecisionValue (id_type new_id) : id(new_id) {}
	
	id_type id;
	
	friend std::size_t hash_value(ESDecisionValue const& b);
	friend ostream& operator<< (ostream& stream, const ESDecisionValue v);
	friend class ESDecisionValueMap;
	friend class ::ESDecisionTreeSuite;
};
std::size_t hash_value(ESDecisionValue const& b);
inline ostream& operator<< (ostream& stream, const ESDecisionValue v){
    return (stream << v.id);
}

/// Value assignment map for ESDecisionValue (manager class).
struct ESDecisionValueMap {
    ESDecisionValueMap () :
	next_bit(0)
    {}
    
    /** Set up a new set of decision values, or confirm they match an existing
     * set (if decision was already entered, and values don't match those
     * associated with the existing decision, a xml_scenario_error is thrown).
     * 
     * @returns The mask covering them all. */
    ESDecisionValue add_decision_values (const string& decision, const std::vector< string > values);
    
    /** Assign from decision and value. add_decision_values must have been
     * called first. */
    ESDecisionValue get (const string& decision, const string& value) const;
    
    typedef map< string, ESDecisionValue > value_map_t;
    /** Get a pair of the following, for decision:
     * first: a mask covering all decision's outputs
     * second: a map of value names to ESDecisionValue objects
     *
     * @throws invalid_argument when decision is not found */
    const pair< ESDecisionValue, value_map_t > getDecision (const string& decision) const;
    
    class ValueFormatter {
	const ESDecisionValueMap& parent;
	const ESDecisionValue value;
	
    public:
	ValueFormatter( const ESDecisionValueMap& p, const ESDecisionValue v ) : parent(p), value(v) {}
	friend ostream& operator<<( ostream& stream, const ValueFormatter& f );
    };
    /** Formats all decision outcomes indicated by an ESDecisionValue, in the
     * format "decision(value), d2(v2)".
     * 
     * Use: "stream << format( value );"
     *
     * This is for error reporting in exceptional situations, and therefore
     * doesn't need to be fast.
     * 
     * This must be a member of ESDecisionValueMap and not ESDecisionValue
     * since ESDecisionValue doesn't know what the codes mean. */
    inline ValueFormatter format( const ESDecisionValue v ) const{
	return ValueFormatter( *this, v );
    }
    /// Implementation of format( ESDecisionValue )
    void format( const ESDecisionValue v, ostream& stream ) const;
    
    private:
	ESDecisionValueMap (const ESDecisionValueMap&) {assert(false);}	// disable copying
	typedef ESDecisionValue::id_type id_type;
	
	// Map of decision to ( pair ( mask, map of value to id ) )
	typedef map< string, pair< ESDecisionValue, value_map_t > > id_map_type;
	id_map_type id_map;
	id_type next_bit;
};
inline ostream& operator<<( ostream& stream, const ESDecisionValueMap::ValueFormatter& f ){
    f.parent.format( f.value, stream );
    return stream;
}

struct ESHostData {
    ESHostData (double aY, WithinHostModel& wH, Pathogenesis::State pS) :
        ageYears(aY), withinHost(wH), pgState(pS) {}
    double ageYears;
    WithinHostModel& withinHost;
    Pathogenesis::State pgState;
};


/** Representation of one decision, random or deterministic (deterministic
 * decisions are hard-coded).
 *
 * Implementations (extending this base) are in the cpp file since they needn't
 * be shared.
 *****************************************************************************/
class ESDecisionTree {
    public:
        /// Run decision tree, with input filtered by mask.
        virtual ESDecisionValue determine (const ESDecisionValue input, const ESHostData& hostData) const =0;
        
	// Note: for some cases we could use ESDecisionName instead of string
	// (for speed), but error messages would be bad. Only slows set-up.
	string decision;	// name of decision
        vector<string> depends;      // other decisions this depends upon
        ESDecisionValue mask;      // mask covering all outputs
        vector<ESDecisionValue> values;    // ids associated with each possible output
        
    protected:
	// Sets mask and values, given the decision's name and a set of values
	void setValues (ESDecisionValueMap& dvMap, const std::vector< string >& valueList);
};

class ESDecisionRandom : public ESDecisionTree {
    public:
	ESDecisionRandom (ESDecisionValueMap& dvMap, const ::scnXml::HSESDecision& xmlDc);
	virtual ESDecisionValue determine (const ESDecisionValue input, const ESHostData& hostData) const;
	
    private:
	// A map from depended decision values (represented as an or-d list of one
	// value (or 0) from each dependency) to a list of cumulative probabilities.
	// Indecies in this list map to the same index in values; last entry must be 1.0.
	//NOTE: be interesting to compare performance with std::map
	typedef unordered_map<ESDecisionValue,vector<double> > map_cum_p_t;
	map_cum_p_t map_cum_p;
	friend struct DR_processor;
};

class ESDecisionUC2Test : public ESDecisionTree {
    public:
	ESDecisionUC2Test (ESDecisionValueMap& dvMap);
	virtual ESDecisionValue determine (const ESDecisionValue, const ESHostData& hostData) const;
};

class ESDecisionAge5Test : public ESDecisionTree {
    public:
	ESDecisionAge5Test (ESDecisionValueMap& dvMap);
	virtual ESDecisionValue determine (const ESDecisionValue input, const ESHostData& hostData) const;
};
class ESDecisionParasiteTest : public ESDecisionTree {
    public:
	ESDecisionParasiteTest (ESDecisionValueMap& dvMap);
	virtual ESDecisionValue determine (const ESDecisionValue input, const ESHostData& hostData) const;
    private:
	// These shouldn't be static: they correspond to the passed-in dvMap and can
	// only be initialized after values are created there.
	// There should only be one ESDecisionParasiteTest object anyway.
	ESDecisionValue test_none, test_microscopy, test_RDT;
};

} }
#endif