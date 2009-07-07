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

/******************************************************************************
 * This is a small program to translate VC's case management trees into a
 * more compact list of outcomes (plus whatever data is needed from the tree).
 * 
 * It's written in D because I (DH) find this the easiest way to program this.
 * See http://www.dsource.org/projects/tango/wiki/TopicInstallTango for setup
 * instructions. I use dsss for building; run: dsss build
 * 
 * API reference:
 * http://www.dsource.org/projects/tango/docs/current/
 *****************************************************************************/
module CaseManagementTree;

import tango.io.device.File;
import tango.io.Stdout;
import tango.util.log.Log;
import tango.util.log.Config;
import tango.text.xml.Document;
import tango.text.xml.DocPrinter;
import tango.math.Math;
import tango.math.IEEE;
import Util=tango.text.Util;
import Float=tango.text.convert.Float;
import Int=tango.text.convert.Integer;

import Clinical.DecisionEnums;

private Logger log;
static this() {
    log = Log.lookup ("CaseManagementTree");
}
/+
/** Possiblity to convert to a smaller representation and check for data-entry errors.
 *
 * Possibly store in files:
 * List of leaf nodes, with attributes like testing method, hospital care or not, delay, treatment plan (small list)
 * Replace decisions list with treatments list, which allow various drug prescription methods or clearing all parasites.
 */

struct Branch {
    int digit;		// digit of decision id adjusted
    double[] nodeProbs;	// probabilities of each outcome; sum(nodeProbs) must add to 1.0
}
struct BranchSet {
    /** Validates this conforms to other branches contained, and adds in to
     * set if probabilities are different from existing members. Returns a
     * reference to the new or existing branch (with same probabilities). */
    Branch* include (Branch b, char[][] bNodeNames, char[] name) {
	with (branches[0]) {	// a branch must exist; we only need to check against one
	    // check format is the same
	    bool valid = true;
	    if (digit != b.digit) valid = false;
	    if (nodeNames.length != bNodeNames.length) valid = false;
	    // extra node vals/IDs are allowed, so only check those present in both
	    for (int i = min(nodeNames.length, bNodeNames.length); i >= 0; --i)
		if (nodeNames[i] != bNodeNames[i]) valid = false;
	    if (!valid)
		throw new Exception ("Branch "~name~" has differing formats");
	}
	foreach (ref existing; branches) {
	    if (existing.nodeProbs.length != b.nodeProbs.length) continue;
	    for (int i = b.nodeProbs.length; i >= 0; --i)
		if (feqrel (existing.nodeProbs[i], b.nodeProbs[i]) < real.mant_dig-1) continue;
	    return &existing;
	}
	branches ~= b;
	return &branches[$-1];
    }
    
    Branch[] branches;
    char[][] nodeNames;	// string value of all outcomes (e.g. "good")
}
BranchSet[char[]] branches;	// set of branch sets, with key branch name (e.g. "quality")
+/

/// This is for filtering out the information we want to pass to the simulator
const DecisionEnums MASK = DecisionEnums.TEST_RDT | DecisionEnums.DRUG_MASK | DecisionEnums.QUALITY_MASK | DecisionEnums.ADHERENCE_MASK | DecisionEnums.MANAGEMENT_MASK | DecisionEnums.TSDELAY_MASK;

int main(char[][] args) {
    if (args.length < 2 || args.length > 3) {
	log.fatal ("Usage: {} infile.xml [outfile.xml]", args[0]);
	return 1;
    }
    char[] inFile = args[1];
    log.info ("Reading file {}", inFile);
    
    auto inDoc = new Document!(char), outDoc = new Document!(char);
    inDoc.parse (cast(char[])File.get (inFile));
    outDoc.header;
    
    auto root = inDoc.elements;
    assert (root.name == "agedependentDecisionTrees");
    auto caseManagements = outDoc.tree.element (null, "caseManagements");
    caseManagements.attribute (null, "name", "NoName");
    
    foreach (decisionTrees; root.query["decisionTrees"].dup) {
	auto caseManagement = caseManagements.element (null, "caseManagement");
	caseManagement.attribute (null, "maxAgeYrs", decisionTrees.query.attribute("maxAgeYrs").nodes[0].value);
	foreach (tree; decisionTrees.query["tree"].dup) {
	    auto caseBranch = caseManagement.element (null, tree.query.attribute("entryPoint").nodes[0].value);
	    
	    double[uint] leaves;	// end points of branches
	    
	    // recurse on "branch"
	    uint endPoints = recurseLeaves (tree, leaves);
	    
	    double prob = 0.0;
	    foreach (p; leaves)
		prob += p;
	    if (feqrel (prob, 1.0) < double.mant_dig-2)	// test prop approx. eq. 1.0
		log.warn ("{} leaves from {} end-points (total probability: {})", leaves.length, endPoints, prob);
	    
	    foreach (id,p; leaves) {
		caseBranch.element (null, "endPoint")
		    .attribute (null, "decision", Int.toString (id, "x#"))
		    .attribute (null, "p", Float.toString (p, double.dig));
	    }
	}
	
	/++ This was intended to be code to read decision stuff - but we
	 + already have enough data to generate this!
	 + Don't even generate it here; do it in the simulator from the id.
	---
	struct DecData {
	    char[] abbrev;
	    enum { POOR, SUFFICIENT, GOOD } qty;
	    int time;
	}
	DecData[int] inDecData;
	auto inDecisions = decisionTrees.query["decisions"].nodes[0];
	foreach (decision; inDecisions.query["decision"].dup) {
	    int id = Int.toInt (decision.query.attribute("id"));
	    assert ((id in inDecData) is null);	// don't want id listed twice
	    auto medicates = decision.query["medicate"];
	    assert (mecicates.count == 1, "code intended to translate from type, simple quantity and time to dosage details");
	    auto med = medicates.nodes[0];
	    DecData dd;
	    dd.abbrev = med.query.attribute("name").nodes[0].value;
	    char[] qty = med.query.attribute("qty").nodes[0].value;
	    if (qty == "good")			dd.qty = GOOD;
	    else if (qty == "sufficient")	dd.qty = SUFFICIENT;
	    else if (qty == "poor")		dd.qty = POOR;
	    else throw new Exception("Expected qty to be good/sufficient/poor");
	    dd.time = med.query.attribute("time").nodes[0].value;
	    inDecData[id] = 
	}
	---
	+/
    }
    
    auto printer = new DocPrinter!(char);
    printer.indent (2);
    if (args.length == 3) {
	File.set (args[2], printer(outDoc));
	log.info ("Written file {}", args[2]);
    } else
	Stdout (printer(outDoc));
    
    // If you want to reformat whitespace in the decision tree file:
    // (compile with version=reformatInput)
    version (reformatInput) {
	printer.cache(false);
	inDoc.header;
	File.set (inFile, printer(inDoc));
    }
    
    return 0;
}

uint recurseLeaves (Document!(char).Node node, ref double[uint] leaves, DecisionEnums inFlags = DecisionEnums.NONE, double inProb = 1.0) {
    auto query = node.query["branch"];
    if (query.count) {
	uint endPoints = 0;
	double cumProb = 0.0;
	foreach (child; query.dup) {
	    DecisionEnums flags = inFlags;
	    char[] decis = child.query.attribute("decision").nodes[0].value;
	    double prob = Float.toFloat(child.query.attribute("p").nodes[0].value);
	    cumProb += prob;
	    prob *= inProb;
	    char[] end;
	    if (prependedBy (decis, "tested:", end)) {
		if (end == "no") {
		} else if (end == "microscopy") {
		    flags |= DecisionEnums.TEST_MICROSCOPY;
		} else if (end == "RDT") {
		    flags |= DecisionEnums.TEST_RDT;
		} else
		    Stdout.format ("unknown - {}", decis).newline;
	    } else if (prependedBy (decis, "result:", end)) {
		if (end == "true positive") {
		    flags |= DecisionEnums.RESULT_TRUE_POS;
		} else if (end == "false negative") {
		    flags |= DecisionEnums.RESULT_FALSE_NEG;
		} else if (end == "false positive") {
		    flags |= DecisionEnums.RESULT_FALSE_POS;
		} else if (end == "true negative") {
		    flags |= DecisionEnums.RESULT_TRUE_NEG;
		} else if (end == "positive") {
		    flags |= DecisionEnums.RESULT_FALSE;
		} else if (end == "negative") {
		    flags |= DecisionEnums.RESULT_TRUE;
		} else
		    Stdout.format ("unknown - {}", decis).newline;
	    } else if (prependedBy (decis, "drug:", end)) {
		if (end == "SP") {
		    flags |= DecisionEnums.DRUG_SP;
		} else if (end == "AL") {
		    flags |= DecisionEnums.DRUG_AL;
		} else if (end == "no antimalarial") {
		    flags |= DecisionEnums.DRUG_NO_AM;
		} else
		    Stdout.format ("unknown - {}", decis).newline;
	    } else if (prependedBy (decis, "quality:", end)) {
		if (end == "good") {
		    flags |= DecisionEnums.QUALITY_GOOD;
		} else if (end == "bad") {
		    flags |= DecisionEnums.QUALITY_BAD;
		} else
		    Stdout.format ("unknown - {}", decis).newline;
	    } else if (prependedBy (decis, "adherence:", end)) {
		if (end == "good") {
		    flags |= DecisionEnums.ADHERENCE_GOOD;
		} else if (end == "bad") {
		    flags |= DecisionEnums.ADHERENCE_BAD;
		} else
		    Stdout.format ("unknown - {}", decis).newline;
	    } else if (prependedBy (decis, "time:", end)) {
		int val = int.min;
		if (end.length == 1)
		    val = end[0] - '0';
		if (val >= 0 || val <= 15) {
		    flags |= (val << DecisionEnums.TSDELAY_LSHIFT);
		} else
		    Stdout.format ("unknown - {}", decis).newline;
	    } else if (prependedBy (decis, "treatment:", end)) {
		if (end == "no antimalarial") {
		    flags |= DecisionEnums.TREATMENT_NO_AM;
		} else if (end == "pre-referral only") {
		    flags |= DecisionEnums.TREATMENT_PREREF;
		} else {
		    flags |= DecisionEnums.DRUG_QN;
		    flags |= DecisionEnums.ADHERENCE_GOOD;
		    if (end == "parenteral") {
			flags |= DecisionEnums.TREATMENT_PARENTAL;
		    } else if (end == "delayed referral") {
			flags |= DecisionEnums.TREATMENT_DELREF | (1 << DecisionEnums.TSDELAY_LSHIFT);
		    } else if (end == "pre-referral and immediate referral") {
			flags |= DecisionEnums.TREATMENT_PREREF_IMMREF;
		    } else if (end == "pre-referral and delayed referral") {
			flags |= DecisionEnums.TREATMENT_PREREF_DELREF | (1 << DecisionEnums.TSDELAY_LSHIFT);
		    } else
			log.error ("unknown: {}", decis);
		}
	    } else if (prependedBy (decis, "management:", end)) {
		if (end == "good") {
		    flags |= DecisionEnums.MANAGEMENT_GOOD;
		} else if (end == "bad") {
		    flags |= DecisionEnums.MANAGEMENT_BAD;
		} else
		    log.error ("unknown: {}", decis);
	    } else {
		log.error ("unknown: {}", decis);
		continue;
	    }
	    endPoints += recurseLeaves (child, leaves, flags, prob);
	}
	if (feqrel (cumProb, 1.0) < double.mant_dig-2) {
	    log.warn ("probabilities of branch [{}] add up to {}",
		node.name == "branch" ?
		    node.query.attribute("decision").nodes[0].value :
		    node.query.attribute("entryPoint").nodes[0].value,
		cumProb);
	}
	return endPoints;
    } else {
	int maskedFlags = MASK & inFlags;
	// leafNode isn't actually required, but for now leave as a placeholder
	query = node.query["leafNode"];
	assert (query.count == 1);
	
	double* pLeaf = maskedFlags in leaves;
	if (pLeaf) {
	    *pLeaf += inProb;
	} else
	    leaves[maskedFlags] = inProb;
	
	return 1;
    }
}

bool prependedBy (char[] source, char[] match, out char[] end) {
    if (source.length >= match.length &&
	source[0..match.length] == match) {
	end = Util.triml(source[match.length..$]);
	return true;
    }
    return false;
}

/** Take val & mask and right-shift by the number of binary zeros less-
 * significant than the first one in mask.
 * 
 * E.g. filterRShift!(0xF0) (0x40) == 0x4. */
uint filterRShift(DecisionEnums mask, uint rShift = 0) (uint val) {
    static assert (rShift < uint.sizeof * 8);
    static if ((mask >> rShift) & 1) {
	return (val & mask) >> rShift;
    } else {
	return filterRShift!(mask, rShift+1) (val);
    }
}
