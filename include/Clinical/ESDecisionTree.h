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
    inline void operator|= (const ESDecisionValue that) {
	id |= that.id;
    }
    private:
	// private constructor, only for use by internal operations
	ESDecisionValue (uint64_t new_id) : id(new_id) {}
	uint64_t id;
	friend std::size_t hash_value(ESDecisionValue const& b);
	friend class ESDecisionValueMap;
};
std::size_t hash_value(ESDecisionValue const& b);

/// Value assignment map for ESDecisionValue (manager class).
struct ESDecisionValueMap {
    ESDecisionValueMap () :
	next_bit(0)
    {}
    
    /** Set up a new set of decision values, returing the mask covering them all. */
    ESDecisionValue add_decision_values (const string& decision, const std::vector< string > values);
    
    /** Assign from decision and value. add_decision_values must have been
     * called first. */
    ESDecisionValue get (const string& decision, const string& value);
    
    private:
	map<string,uint64_t> id_map;
	uint64_t next_bit;
};

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
        virtual ESDecisionValue determine (const ESDecisionValue input, ESHostData& hostData) const =0;
        
	// Note: for some cases we could use ESDecisionName instead of string
	// (for speed), but error messages would be bad. Only slows set-up.
	string decision;	// name of decision
        vector<string> depends;      // other decisions this depends upon
        ESDecisionValue mask;      // mask covering all outputs
        vector<ESDecisionValue> values;    // ids associated with each possible output
        
    protected:
	// Sets mask and values, given the decision's name and a set of values
	void setValues (ESDecisionValueMap dvMap, const std::vector< string >& valueList);
};

class ESDecisionRandom : public ESDecisionTree {
    public:
	ESDecisionRandom (ESDecisionValueMap dvMap, const ::scnXml::Decision& xmlDc);
	virtual ESDecisionValue determine (const ESDecisionValue input, ESHostData& hostData) const;
	
    private:
	// a map of depended decision values to
	// cumulative probabilities associated with values
	//NOTE: be interesting to compare performance with std::map
	unordered_map<ESDecisionValue,vector<double> > set_cum_p;
};

class ESDecisionUC2Test : public ESDecisionTree {
    public:
	ESDecisionUC2Test (ESDecisionValueMap dvMap);
	virtual ESDecisionValue determine (const ESDecisionValue, ESHostData& hostData) const;
};

class ESDecisionAge5Test : public ESDecisionTree {
    public:
	ESDecisionAge5Test (ESDecisionValueMap dvMap);
	virtual ESDecisionValue determine (const ESDecisionValue input, ESHostData& hostData) const;
};
class ESDecisionParasiteTest : public ESDecisionTree {
    public:
	ESDecisionParasiteTest (ESDecisionValueMap dvMap);
	virtual ESDecisionValue determine (const ESDecisionValue input, ESHostData& hostData) const;
    private:
	// These shouldn't be static: they correspond to the passed-in dvMap and can
	// only be initialized after values are created there.
	// There should only be one ESDecisionParasiteTest object anyway.
	ESDecisionValue test_none, test_microscopy, test_RDT;
};

} }
#endif