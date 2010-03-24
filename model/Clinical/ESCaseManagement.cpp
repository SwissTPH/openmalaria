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

#include "Clinical/ESCaseManagement.h"
#include "inputData.h"
#include "util/random.h"
#include "util/errors.hpp"

#include <boost/format.hpp>
#include <boost/static_assert.hpp>
#include <boost/unordered_map.hpp>
#include <boost/tokenizer.hpp>

namespace OM { namespace Clinical {
    using namespace OM::util;

// -----  Types  -----

void ESDecisionTree::setValues (const char* decision, vector<const char*> valueList) {
    mask = ESDecisionValue::add_decision_values (decision, valueList);
    values.resize (valueList.size());
    size_t i = 0;
    BOOST_FOREACH ( const char* value, valueList ) {
	values[i].assign (decision, value);
    }
}
/*class ESDecisionDeterministic : public ESDecisionTree {
public:
    
    virtual ESDecisionValue determine (ESDecisionValue input, ESHostData& hostData) {
   unordered_map<ESDecisionValue,ESDecisionValue>::const_iterator it = set.find (input);
   if (it == set.end())
   return 0;
   return it->second;
}

    private:
	unordered_map<ESDecisionValue,ESDecisionValue> set;
};*/
class ESDecisionRandom : public ESDecisionTree {
    public:
	ESDecisionRandom (const ::scnXml::Decision& xmlDc) {
	    // Set depends, values:
	    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
	    boost::char_separator<char> separator (",");
	    
	    tokenizer depends_tokens (xmlDc.getDepends(), separator);
	    BOOST_FOREACH ( const string& str, depends_tokens ) {
		//FIXME: check names confirm to expected format
		depends.push_back (str.c_str());	// cast to ESDecisionName
	    }
	    
	    vector<const char*> valueList;
	    tokenizer values_tokens (xmlDc.getValues(), separator);
	    BOOST_FOREACH ( const string& str, values_tokens ) {
		valueList.push_back (str.c_str());
	    }
	    setValues (xmlDc.getName().c_str(), valueList);
	    
	    //TODO: parse tree
	}
	
	virtual ESDecisionValue determine (ESDecisionValue input, ESHostData& hostData) {
	    unordered_map<ESDecisionValue,vector<double> >::const_iterator it = set_cum_p.find (input);
	    if (it == set_cum_p.end())
		return ESDecisionValue();	// no decision
	    double sample = random::uniform_01 ();
	    size_t i = 0;
	    while (it->second[i] < sample)
		++i;
	    return values[i];
	}
	
    private:
	// a map of depended decision values to
	// cumulative probabilities associated with values
	//NOTE: be interesting to compare performance with std::map
	unordered_map<ESDecisionValue,vector<double> > set_cum_p;
};
class ESDecisionAge5Test : public ESDecisionTree {
    public:
	ESDecisionAge5Test () {
	    vector<const char*> valueList (2, "under5");
	    valueList[1] = "over5";
	    setValues ("age5Test", valueList);
	}
	virtual ESDecisionValue determine (ESDecisionValue input, ESHostData& hostData) {
	    if (hostData.ageYears >= 5.0)
		return values[1];
	    else
		return values[0];
	}
};
class ESDecisionParasiteTest : public ESDecisionTree {
    public:
	ESDecisionParasiteTest () {
	    depends.resize (1, "test");
	    //TODO: assert "test" has outputs as below
	    
	    vector<const char*> valueList (2, "negative");
	    valueList[1] = "positive";
	    setValues ("result", valueList);
	}
	
	virtual ESDecisionValue determine (ESDecisionValue input, ESHostData& hostData) {
	static const ESDecisionValue TEST_MICROSCOPY("test", "microscopy"),
		TEST_RDT("test","RDT"),
		TEST_NONE("test","none");
	    if (input == TEST_MICROSCOPY) {
		// if (hostData.withinHost.getDensity () > ...)
		//FIXME: or values[0]
		return values[1];
	    } else if (input == TEST_RDT) {
		//FIXME: or values[0]
		return values[0];
	    } else
		return ESDecisionValue();	// 0, no decision
	}
};


// ESCaseManagement::TreeType ESCaseManagement::cmTree;
// cmid ESCaseManagement::cmMask;

ESDecisionMap ESCaseManagement::UC, ESCaseManagement::severe;
CaseTreatment* ESCaseManagement::mdaDoses;

// -----  Static  -----

void ESCaseManagement::init () {
    // Assume EventScheduler data was checked present:
    const scnXml::HSEventScheduler& xmlESCM = InputData().getHealthSystem().getEventScheduler().get();
    UC.initialize (xmlESCM.getUC (), false);
    severe.initialize (xmlESCM.getSevere (), true);
    
    // MDA Intervention data
    const scnXml::Interventions::MDADescriptionOptional mdaDesc = InputData().getInterventions().getMDADescription();
    if (mdaDesc.present()) {
	mdaDoses = new CaseTreatment (/*FIXME mdaDesc.get().getMedicate()*/);
    } else
	mdaDoses = NULL;
}

void ESDecisionMap::initialize (const ::scnXml::HSESCMDecisions& xmlDcs, bool complicated) {
    uint32_t decision_num = xmlDcs.getDecision().size() + 1;
    if (!complicated)
	decision_num += 1;
    decisions.resize (decision_num);
    
    decision_num = 0;
    decisions[decision_num++] = new ESDecisionAge5Test;
    if (!complicated) {
	decisions[decision_num++] = new ESDecisionParasiteTest;
    }
    
    BOOST_FOREACH ( const ::scnXml::Decision& xmlDc, xmlDcs.getDecision() ) {
	decisions[decision_num++] = new ESDecisionRandom (xmlDc);
    }
}
ESDecisionMap::~ESDecisionMap () {
    BOOST_FOREACH ( ESDecisionTree* d, decisions ) {
	delete d;
    }
}


void ESCaseManagement::massDrugAdministration(list<MedicateData>& medicateQueue) {
    if (mdaDoses == NULL)
	throw util::xml_scenario_error ("MDA intervention without description");
    else
	mdaDoses->apply(medicateQueue);
}

void ESDecisionName::operator= (const char* name) {
    map<string,uint32_t>::const_iterator it = id_map.find (name);
    if (it == id_map.end()) {
        id = next_free;
        next_free++;
        id_map[name] = id;
    } else {
        id = it->second;
    }
}
map<string,uint32_t> ESDecisionName::id_map;
uint32_t ESDecisionName::next_free = 1;

ESDecisionValue ESDecisionValue::add_decision_values (const char* decision, vector<const char*> values) {
    // got length l = values.size() + 1 (default, "no outcome"); want minimal n such that: 2^n >= l
    // that is, n >= log_2 (l)
    // so n = ceil (log_2 (l))
    uint32_t n_bits = std::ceil (log (values.size() + 1) / log(2.0));
    if (n_bits+next_bit>=(sizeof(next_bit)*8))
	throw runtime_error ("ESDecisionValue design: insufficient bits");
    uint64_t next=(1<<next_bit), step;
    step=next;
    BOOST_FOREACH ( const char* value, values ) {
	string name = (boost::format("%1%(%2%)") %decision %value).str();
	bool inserted = id_map.insert (pair<string,uint64_t>(name, next)).second;
	if (!inserted)
	    throw runtime_error ("ESDecisionValue: value doesn't have a unique name");
	next += step;
    }
    assert (next <= (1<<(n_bits+next_bit)));
    
    ESDecisionValue mask;
    for (size_t i = 0; i < n_bits; ++i) {
	mask.id |= (1<<next_bit);
	++next_bit;
    }
    assert ((next-step) & mask.id == 0);
    assert (next & mask.id != 0);
    return mask;
}
void ESDecisionValue::assign (const char* decision, const char* value) {
    string name = (boost::format("%1%(%2%)") %decision %value).str();
    map<string,uint64_t>::const_iterator it = id_map.find (name);
    if (it == id_map.end()) {
	throw runtime_error ("ESDecisionValue: name used before an entry was created for it");
    } else {
	id = it->second;
    }
}
std::size_t hash_value(ESDecisionValue const& b) {
    boost::hash<int> hasher;
    return hasher(b.id);
}
map<string,uint64_t> ESDecisionValue::id_map;
uint64_t ESDecisionValue::next_bit = 0;

CaseTreatment* ESDecisionMap::determine (OM::Clinical::ESHostData& hostData) {
    ESDecisionValue id;
    //TODO: add pgState decision
    //TODO
}


void ESCaseManagement::execute (list<MedicateData>& medicateQueue, Pathogenesis::State pgState, WithinHost::WithinHostModel& withinHostModel, double ageYears, SurveyAgeGroup ageGroup) {
    ESDecisionMap* map;
    assert (pgState & Pathogenesis::SICK);
    if (pgState & Pathogenesis::COMPLICATED)
        map = &severe;
    else
        map = &UC;
    
    ESHostData hostData (ageYears, withinHostModel, pgState);
    
    CaseTreatment* treatment = map->determine (hostData);
    
    // We always remove any queued medications.
    medicateQueue.clear();
    treatment->apply (medicateQueue);
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