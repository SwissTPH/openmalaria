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

#include <boost/static_assert.hpp>
#include <boost/unordered_map.hpp>

namespace OM { namespace Clinical {
    using namespace OM::util;

// -----  Types  -----

/*class ESDecisionDeterministic : public ESDecisionTree {
public:
    
    virtual ESDecisionID determine (ESDecisionID input, ESHostData& hostData) {
   unordered_map<ESDecisionID,ESDecisionID>::const_iterator it = set.find (input);
   if (it == set.end())
   return 0;
   return it->second;
}

    private:
	unordered_map<ESDecisionID,ESDecisionID> set;
};*/
class ESDecisionRandom : public ESDecisionTree {
    public:
	ESDecisionRandom (const ::scnXml::Decision& xmlDc) {
	    //TODO: set depends, values
	    //TODO: parse tree
	}
	
	virtual ESDecisionID determine (ESDecisionID input, ESHostData& hostData) {
	    unordered_map<ESDecisionID,vector<double> >::const_iterator it = set_cum_p.find (input);
	    if (it == set_cum_p.end())
		return 0;
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
	unordered_map<ESDecisionID,vector<double> > set_cum_p;
	vector<ESDecisionTree*> decisions;
};
class ESDecisionAge5Test : public ESDecisionTree {
    public:
	ESDecisionAge5Test () {
	    mask = Decision::AGE_OVER5 | Decision::AGE_UNDER5;
	    values.resize (2, Decision::AGE_OVER5);
	    values[1] = Decision::AGE_UNDER5;
	}
	virtual ESDecisionID determine (ESDecisionID input, ESHostData& hostData) {
	    if (hostData.ageYears >= 5.0)
		return Decision::AGE_OVER5;
	    else
		return Decision::AGE_UNDER5;
	}
};
class ESDecisionParasiteTest : public ESDecisionTree {
    public:
	ESDecisionParasiteTest () {
	    depends.resize (1, "tested");
	    mask = Decision::RESULT_MASK;
	    values.resize (2, Decision::RESULT_POSITIVE);
	    values[1] = Decision::RESULT_NEGATIVE;
	}
	
	virtual ESDecisionID determine (ESDecisionID input, ESHostData& hostData) {
	    if (input == Decision::TEST_MICROSCOPY) {
		// if (hostData.withinHost.getDensity () > ...)
		//FIXME: or RESULT_NEGATIVE
		return Decision::RESULT_POSITIVE;
	    } else if (input == Decision::TEST_RDT) {
		//FIXME: or RESULT_NEGATIVE
		return Decision::RESULT_POSITIVE;
	    } else
		return 0;
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

CaseTreatment* ESDecisionMap::determine (Pathogenesis::State pgState, ESHostData& hostData) {
    BOOST_STATIC_ASSERT ( int(Pathogenesis::MORBIDITY_MASK) == int(Decision::MORBIDITY_MASK) );
    ESDecisionID id = pgState & Decision::MORBIDITY_MASK;
    //TODO
}


void ESCaseManagement::execute (list<MedicateData>& medicateQueue, Pathogenesis::State pgState, WithinHost::WithinHostModel& withinHostModel, double ageYears, SurveyAgeGroup ageGroup) {
    ESDecisionMap* map;
    assert (pgState & Pathogenesis::SICK);
    if (pgState & Pathogenesis::COMPLICATED)
        map = &severe;
    else
        map = &UC;
    
    ESHostData hostData (ageYears, withinHostModel);
    
    CaseTreatment* treatment = map->determine (pgState, hostData);
    
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