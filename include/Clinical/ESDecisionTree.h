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

#ifndef Hmod_ESDecisionTree
#define Hmod_ESDecisionTree

#include "Global.h"
#include "Clinical/ESDecisionValue.h"
#include "Pathogenesis/State.h"
#include "WithinHost/WithinHostModel.h"

#include <map>
#include <boost/unordered_map.hpp>

namespace scnXml{
    class HSESDecision;
}
namespace OM { namespace Clinical {
    using WithinHost::WithinHostModel;
    using boost::unordered_map;

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
	/// Create a user-configured decision from xmlDc.
	static ESDecisionTree* create (ESDecisionValueMap& dvm, const ::scnXml::HSESDecision& xmlDc);
	
        /** Run decision tree. */
	inline ESDecisionValue determine (const ESDecisionValue input, const ESHostData& hostData) const{
	    return determineImpl( input & mask, hostData );
	}
	
	string decision;	// name of decision
	vector<string> depends;      // other decisions this depends upon
	
    protected:
        virtual ESDecisionValue determineImpl (const ESDecisionValue input, const ESHostData& hostData) const =0;
        ESDecisionValue mask;      // mask covering all dependencies' values
};

class ESDecisionAge : public ESDecisionTree {
    public:
	ESDecisionAge (ESDecisionValueMap& dvMap, const ::scnXml::HSESDecision& xmlDc);
	virtual ESDecisionValue determineImpl (const ESDecisionValue input, const ESHostData& hostData) const;
	
    private:
	// A map from age-group upper-bounds to output values. The first lower-
	// bound we assume to be less than or equal to any input value, and we
	// assume no input can be greater than the last upper-bound (which
	// should be inf).
	map<double, ESDecisionValue> age_upper_bounds;
	
	friend struct DA_processor;
};

class ESDecisionRandom : public ESDecisionTree {
    public:
	/** Set up.
	 *
	 * @param dvMap Decision-value map to add decision-outcomes into.
	 * @param xmlDc XML element describing tree
	 * @param dependsInput Prerequisites of this decision. (Also described
	 * within xmlDc; passed to avoid re-parsing.) */
	ESDecisionRandom (ESDecisionValueMap& dvMap, const ::scnXml::HSESDecision& xmlDc, const vector<string>& dependsInput);
	virtual ESDecisionValue determineImpl (const ESDecisionValue input, const ESHostData& hostData) const;
	
    private:
	// A map from depended decision values (represented as an or-d list of one
	// value (or 0) from each dependency) to a list of cumulative probabilities.
	// Indecies in this list map to the same index in values; last entry must be 1.0.
	//NOTE: be interesting to compare performance between boost::unordered_map and std::map
	typedef unordered_map<ESDecisionValue,vector<double> > map_cum_p_t;
	map_cum_p_t map_cum_p;
	
	vector<ESDecisionValue> values;    // ids associated with each possible output
	
	friend struct DR_processor;
};

class ESDecisionUC2Test : public ESDecisionTree {
    public:
	ESDecisionUC2Test (ESDecisionValueMap& dvMap);
	virtual ESDecisionValue determineImpl (const ESDecisionValue, const ESHostData& hostData) const;
    private:
	ESDecisionValue UC1, UC2;
};

class ESDecisionParasiteTest : public ESDecisionTree {
    public:
	ESDecisionParasiteTest (ESDecisionValueMap& dvMap);
	virtual ESDecisionValue determineImpl (const ESDecisionValue input, const ESHostData& hostData) const;
    private:
	ESDecisionValue none, positive, negative;
	// These shouldn't be static: they correspond to the passed-in dvMap and can
	// only be initialized after values are created there.
	// There should only be one ESDecisionParasiteTest object anyway.
	ESDecisionValue test_none, test_microscopy, test_RDT;
};

} }
#endif