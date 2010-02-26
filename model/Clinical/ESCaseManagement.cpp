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
#include "util/gsl.h"
#include "util/errors.hpp"

#include <boost/static_assert.hpp>

namespace OM { namespace Clinical {

ESCaseManagement::TreeType ESCaseManagement::cmTree;
cmid ESCaseManagement::cmMask;

CaseTreatment* ESCaseManagement::mdaDoses;

// -----  Static  -----

void ESCaseManagement::init () {
    const scnXml::CaseManagementTree& xmlCM = InputData.getEventScheduler().getCaseManagementTree();
    for (scnXml::CaseManagementTree::CM_pBranchSetConstIterator it = xmlCM.getCM_pBranchSet().begin(); it != xmlCM.getCM_pBranchSet().end(); ++it) {
	cmTree[it->getID ()] = new CMPBranchSet (it->getCM_pBranch());
    }
    for (scnXml::CaseManagementTree::CM_leafConstIterator it = xmlCM.getCM_leaf().begin(); it != xmlCM.getCM_leaf().end(); ++it) {
	cmTree[it->getID ()] = new CMLeaf (CaseTreatment (it->getMedicate()));
    }
    
    cmMask = xmlCM.getMask();
    
    // MDA Intervention data
    const scnXml::Interventions::MDADescriptionOptional mdaDesc = InputData.getInterventions().getMDADescription();
    if (mdaDesc.present()) {
	mdaDoses = new CaseTreatment (mdaDesc.get().getMedicate());
    } else
	mdaDoses = NULL;
}


void ESCaseManagement::massDrugAdministration(list<MedicateData>& medicateQueue) {
    if (mdaDoses == NULL)
	throw util::xml_scenario_error ("MDA intervention without description");
    else
	mdaDoses->apply(medicateQueue, 0);
}


cmid ESCaseManagement::execute (list<MedicateData>& medicateQueue, Pathogenesis::State pgState, WithinHost::WithinHostModel& withinHostModel, double ageYears, SurveyAgeGroup ageGroup) {
#ifndef NDEBUG
    if (! (pgState & Pathogenesis::SICK))
        throw domain_error ("doCaseManagement shouldn't be called if not sick");
#endif

    // We always remove any queued medications.
    medicateQueue.clear();
    
    BOOST_STATIC_ASSERT ( cmid(Pathogenesis::MORBIDITY_MASK) == cmid(Decision::MORBIDITY_MASK) );
    cmid decisionID = pgState & Pathogenesis::MORBIDITY_MASK;;
    if (ageYears > 5.0)
	decisionID |= Decision::AGE_OVER5;
    
    CaseTreatmentPair leaf = traverse (decisionID);	// Traverse the tree.
    leaf.second.apply (medicateQueue, leaf.first);
    return leaf.first;
}

CaseTreatmentPair ESCaseManagement::traverse (cmid id) {
    /* This has branches for severe or uncomplicated cases.
    FIXME: In UCs it needs to branch into UC1 or UC2 depending on
    result of last parasite test, if within health-system-memory. */
    
    // Parasite test:
    //NOTE: currently only one test, but only used in UC trees
    // (can make test depend on severe/UC state)
    if ((id & Decision::RESULT_MASK) == Decision::RESULT_DETERMINE) {
	id = id & (~Decision::RESULT_MASK);	// remove result flags
	
	//FIXME: or RESULT_NEGATIVE
	id |= Decision::RESULT_POSITIVE;
	
	return traverse(id);
    }
    
    TreeType::const_iterator it = cmTree.find (id & cmMask);
    if (it == cmTree.end()) {
	ostringstream msg;
	msg << "No node for id "<<(id&cmMask)<<" (unmasked: "<<id<<")";
	throw util::xml_scenario_error (msg.str());
    }
    return it->second->traverse (id);
}


// -----  CMNode derivatives  -----

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
}

CaseTreatmentPair ESCaseManagement::CMPBranchSet::traverse (cmid id) {
    double randCum = rng::uniform01();
    size_t i = 0;
    while (branches[i].cumP < randCum) {
	++i;
	assert (i < branches.size());
    }
    id |= branches[i].outcome;
    
    return ESCaseManagement::traverse (id);
}

CaseTreatmentPair ESCaseManagement::CMLeaf::traverse (cmid id) {
    return CaseTreatmentPair (id, ct);
}

} }