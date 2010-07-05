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

// Utility code to do with representing decisions.

#ifndef Hmod_ESDecisionValue
#define Hmod_ESDecisionValue

#include "Global.h"

#include <cassert>
#include <boost/tuple/tuple.hpp>

class ESDecisionTreeSuite;

namespace OM { namespace Clinical {
    using boost::tuple;

/** A compressed representation of all decision outcomes.
 * Pass by value (it is 64-bits in size). */
struct ESDecisionValue {
    ESDecisionValue () : id(0) {}
    
    inline bool operator== (const ESDecisionValue that) const {
	return id == that.id;
    }
    inline bool operator!= (const ESDecisionValue that) const {
	return id != that.id;
    }
    inline ESDecisionValue operator& (const ESDecisionValue that) const {
	return ESDecisionValue(id & that.id);
    }
    inline ESDecisionValue operator| (const ESDecisionValue that) const {
	return ESDecisionValue( id | that.id );
    }
    inline void operator|= (const ESDecisionValue that) {
	id |= that.id;
    }
    // all use in a map:
    inline bool operator< (const ESDecisionValue that) const {
	return id < that.id;
    }
    
private:
    typedef boost::uint64_t id_type;
    
    // private constructor, only for use by internal operations
    ESDecisionValue (id_type new_id) : id(new_id) {}
    
    id_type id;
    
    friend std::size_t hash_value(ESDecisionValue const& b);
    friend ostream& operator<< (ostream& stream, const ESDecisionValue v);
    friend struct ESDecisionValueMap;
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
    /** Reset to zero. */
    void clear() {
	id_map.clear();
	next_bit = 0;
    }
    
    /** Set up a new set of decision values, or confirm they match an existing
     * set (if decision was already entered, and values don't match those
     * associated with the existing decision, a xml_scenario_error is thrown).
     * 
     * To allow setting up bit-masks, values are assigned integer values in the
     * order added (as 0, s, 2s, etc., where s is some step size). This only
     * applies when decision wasn't already added; its recommended to test
     * values with assert.
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
    tuple< ESDecisionValue, const value_map_t& > getDecision (const string& decision) const;
    
    /** Similar to getDecision, but just get the mask directly. */
    ESDecisionValue getDecisionMask (const string& decision) const;
    
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

/** Encapsulates a set of ESDecisionValues, allowing piecewise operators on the
 * set. */
struct ESDecisionValueSet {
    /// Construct with just one value: [ val ]
    inline ESDecisionValueSet (ESDecisionValue val) {
	values.push_back( val );
    }
    
    /// Construct from all values of a value_map_t map
    ESDecisionValueSet (const ESDecisionValueMap::value_map_t valMap);
    
    /// This set becomes the set of all elements x | y such that x is in this
    /// and y is in that.
    void operator|= (const ESDecisionValueSet& that);
    
    list<ESDecisionValue> values;
};

} }
#endif
