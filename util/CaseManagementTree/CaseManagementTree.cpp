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

/******************************************************************************
 * This is a small program to translate VC's case management trees into a
 * more compact list of outcomes (plus whatever data is needed from the tree).
 * 
 * This is a C++ conversion of the D program to avoid dependance on D.
 *****************************************************************************/

#include "Clinical/ESDecision.h"
#include "Pathogenesis/State.h"
using namespace OM;
using namespace OM::Clinical;

// Use glib's libxml2 here; this is not the same xml library used by the openmalaria simulator but
// (a) this script won't be compiled a lot so it doesn't matter and (b) I wanted to evaluate the
// usage of libxml2 (or C++ wrapper libxml++).
#include <libxml++/libxml++.h>
using Glib::ustring;

#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <assert.h>

using namespace std;
using namespace boost;

typedef Decision::DecisionEnums d_id_t;

/// This is for filtering out the information we want to pass to the simulator
const d_id_t MASK = d_id_t (Decision::TEST_RDT | Decision::DRUG_MASK | Decision::QUALITY_MASK | Decision::ADHERENCE_MASK | Decision::TREATMENT_MASK | Decision::TSDELAY_MASK);

/// Function to initialize decisionsMap
const map<ustring,map<ustring,d_id_t> > createDecisionsMap () {
    map<ustring,map<ustring,d_id_t> > decisionsMap;
    map<ustring,d_id_t>* currentMap = NULL;
    
    currentMap = &decisionsMap["maxAge"];
    (*currentMap)["5"] = Decision::NONE;
    (*currentMap)["99"] = Decision::AGE_OVER5;
    
    //TODO: severe & second case?
    currentMap = &decisionsMap["case"];
    (*currentMap)["UC1"] = d_id_t (Pathogenesis::SICK);
    (*currentMap)["UC2"] = d_id_t (Pathogenesis::SICK | Pathogenesis::SECOND_CASE);
    (*currentMap)["severe"] = d_id_t (Pathogenesis::SICK | Pathogenesis::SEVERE);
    
    currentMap = &decisionsMap["source"];
    (*currentMap)["hospital"] = Decision::TREATMENT_HOSPITAL;
    
    currentMap = &decisionsMap["tested"];
    (*currentMap)["none"] = Decision::TEST_NONE;
    (*currentMap)["microscopy"] = Decision::TEST_MICROSCOPY;
    (*currentMap)["RDT"] = Decision::TEST_RDT;
    
    currentMap = &decisionsMap["result"];
    (*currentMap)["positive"] = Decision::RESULT_POSITIVE;
    (*currentMap)["negative"] = Decision::RESULT_NEGATIVE;
    // All we're interested in right here is the test outcome — right?
    (*currentMap)["true positive"] = Decision::RESULT_POSITIVE;
    (*currentMap)["false negative"] = Decision::RESULT_NEGATIVE;
    
    currentMap = &decisionsMap["drug"];
    (*currentMap)["no antimalarial"] = Decision::DRUG_NO_AM;
    (*currentMap)["SP"] = Decision::DRUG_SP;
    (*currentMap)["AL"] = Decision::DRUG_AL;
    
    currentMap = &decisionsMap["adherence"];
    (*currentMap)["good"] = Decision::ADHERENCE_FULL;
    (*currentMap)["missed first dose"] = Decision::ADHERENCE_MISSED_FIRST;
    //TODO: not the same as ADHERENCE_MISSED_LAST?
//     (*currentMap)["missed last day"] = Decision::ADHERENCE_MISSED_LAST_DAY;
    //TODO: remove (added for compatibility with old tree)
    (*currentMap)["bad"] = Decision::ADHERENCE_SELECTIVE;
    
    currentMap = &decisionsMap["quality"];
    (*currentMap)["good"] = Decision::QUALITY_GOOD;
    (*currentMap)["bad"] = Decision::QUALITY_BAD;
    
    currentMap = &decisionsMap["time"];
    (*currentMap)["0"] = Decision::NONE;
    (*currentMap)["1"] = d_id_t (1 << Decision::TSDELAY_SHIFT);
    (*currentMap)["2"] = d_id_t (2 << Decision::TSDELAY_SHIFT);
    
    return decisionsMap;
}

/** Map of decisions to maps of value (at decision) to id.
 *
 * Initialized once through a neat trick (function call converts to const). */
const map<ustring,map<ustring,d_id_t> > decisionsMap = createDecisionsMap ();

/** Function to safely get a value from decisionsMap. */
const map<ustring,d_id_t>& decisionsMapGet (const ustring& k) {
    const map<ustring,map<ustring,d_id_t> >::const_iterator it = decisionsMap.find (k);
    if (it == decisionsMap.end()) {
	ostringstream msg;
	msg << "depends \"" << k << "\" unrecognized";
	throw xmlpp::exception (msg.str());
    }
    return it->second;
}


//TODO: build tree
//TODO: we could probably often bubble input-decisions up over random decisions and maybe sometimes
//the other way. If this reduces the number of random trees required it is an advantage.


/// First-stage parsing: produce an in-memory tree. This is only required to
/// deal with id/ref of branches.
class CMRefTreeParser : public xmlpp::SaxParser
{
public:
    struct TreeNode {
	virtual ~TreeNode() {}
    };
    struct TreeBranches;
    struct TreeChoice : public TreeNode {
	TreeChoice () : child(NULL), id(Decision::NONE), prob(numeric_limits<double>::signaling_NaN()) {}
	TreeChoice (const TreeBranches& parent, const Glib::ustring& value) : child(NULL), id((d_id_t)0xFFFFFFFF), prob(numeric_limits<double>::signaling_NaN()) {
	    map<ustring,d_id_t>::const_iterator local_id = parent.id_value_map.find (value);
	    if (local_id == parent.id_value_map.end()) {
		ostringstream msg;
		msg << "unexpected choice value: " << value;
		throw xmlpp::exception (msg.str());
	    }
	    id = local_id->second;
	};
	
	TreeBranches* child;
	d_id_t id;
	double prob;
    };
    struct TreeBranches : public TreeNode {
	TreeBranches (const TreeChoice& parent, const Glib::ustring& depends) :
	id_value_map(decisionsMapGet(depends)),
	local_cum_prob (0.0)
	{}
	
	list<TreeChoice> choices;	// all choices possible for this branch-point
	const map<ustring,d_id_t>& id_value_map;	// map to resolve an id from a value, for this decision
	double local_cum_prob;	// initialized to zero and incermented for each choice; should come to 1.0
    };
    
    TreeChoice root;
    
    
    CMRefTreeParser() : line_num(1) {}
//     virtual ~CMRefTreeParser();
    
protected:
    //overrides:
//     virtual void on_start_document() {}
//     virtual void on_end_document();
    virtual void on_start_element(const Glib::ustring& name, const AttributeList& attributes) {
	try {
	    if (name == "agedependentDecisionTrees") {
		confirm (elt_stack.empty (), "agedependentDecisionTrees should only be first (root) element");
		
		// root.id left at initialization value
		root.prob = 1.0;
		elt_stack.push_back (make_pair(name,&root));
	    } else {
		confirm (!elt_stack.empty (), "expected <agedependentDecisionTrees> as first (root) element");
		
		if (name == "randomBranches" || name == "inputBranches") {
		    ustring attrKey = "ref";
		    AttributeList::const_iterator attribute = std::find_if(attributes.begin(), attributes.end(), AttributeHasName(attrKey));
		    if (attribute != attributes.end()) {	// is a "ref" to another *Branches
			ostringstream msg;
			msg << name << " reference not found: " << attribute->value;
			confirm (refMap.count(attribute->value), msg.str());
			elt_stack.push_back (make_pair(name, refMap[attribute->value]));
		    } else {
			TreeChoice* parent = dynamic_cast<TreeChoice*> (elt_stack.back().second);
			confirm (parent != NULL, "*Branches should only be a child of a choice (or agedependentDecisionTrees) element");
			
			TreeBranches* branches = new TreeBranches (*parent, get_attribute (attributes, "depends", name));
			parent->child = branches;
			elt_stack.push_back (make_pair (name, branches));
			
			attrKey = "id";
			AttributeList::const_iterator attribute = std::find_if(attributes.begin(), attributes.end(), AttributeHasName(attrKey));
			if (attribute != attributes.end()) {
			    refMap[attribute->value] = branches;
			}
		    }
		} else if (name == "choice") {
		    TreeBranches* parent = dynamic_cast<TreeBranches*> (elt_stack.back().second);
		    confirm (parent != NULL, "choice should only be a child of a *Branches element");
		    
		    parent->choices.push_back (TreeChoice (*parent, get_attribute (attributes, "value", name)));
		    TreeChoice* choice = &parent->choices.back();
		    
		    if (elt_stack.back().first == "randomBranches") {
			double p = lexical_cast<double> (get_attribute (attributes, "p", "choice element when inside a randomBranches"));
			parent->local_cum_prob += p;
			choice->prob = p;
		    }
		    
		    elt_stack.push_back (make_pair(name, choice));
		}
	    }
	} catch (xmlpp::exception& e) {
	    cerr << "line " << line_num << ": ";
	    throw;
	}
    }
    virtual void on_end_element(const Glib::ustring& name) {
	try {
	    if (elt_stack.empty() || name != elt_stack.back().first) {
		ostringstream msg;
		msg << "mismatched tags: ";
		if (elt_stack.empty()) {
		    msg << "(none)";
		} else {
		    msg << '<' << elt_stack.back().first << '>';
		}
		msg << " and </" << name << '>';
		throw xmlpp::exception (msg.str());
	    }
	    if (elt_stack.back().first == "randomBranches") {
		TreeBranches* parent = dynamic_cast<TreeBranches*> (elt_stack.back().second);
		confirm (parent != NULL, "element should be a TreeBranches");
		if (!(0.999 < parent->local_cum_prob && parent->local_cum_prob < 1.001)) {
		    ostringstream msg;
		    msg << "probabilities of randomBranches's children should add up to 1.0, not " << parent->local_cum_prob;
		    throw xmlpp::exception (msg.str());
		}
	    } else if (elt_stack.back().first == "choice") {
		TreeChoice* choice = dynamic_cast<TreeChoice*> (elt_stack.back().second);
		confirm (choice != NULL, "element should be a TreeChoice");
		if (choice->child == NULL)
		    cout << format("choice: %1$_#0.9x\tp: %2%") % choice->id % choice->prob <<endl;
	    }
	    elt_stack.pop_back();
	} catch (xmlpp::exception& e) {
	    cerr << "line " << line_num << ": ";
	    throw;
	}
    }
    //NOTE: we try to count new-lines here; doesn't seem to be quite accurate
    virtual void on_characters(const Glib::ustring& characters) {
	for (Glib::ustring::const_iterator it = characters.begin(); it != characters.end(); ++it) {
	    gunichar c = *it;
	    if (c == '\n')
		line_num++;
	    else if (!Glib::Unicode::isspace (c)) {
		cerr << "line " << line_num << ": ";
		throw xmlpp::exception ("Character data not expected");
	    }
	}
    }
    virtual void on_comment(const Glib::ustring& text) {
	for (Glib::ustring::const_iterator it = text.begin(); it != text.end(); ++it) {
	    if (*it == '\n')
		line_num++;
	}
    }
    virtual void on_warning(const Glib::ustring& text) {
	cerr << "line " << line_num << ": ";
	cerr << "Warning: " << text << endl;
    }
    virtual void on_error(const Glib::ustring& text) {
	cerr << "line " << line_num << ": ";
	cerr << "Error: " << text << endl;
    }
    virtual void on_fatal_error(const Glib::ustring& text) {
	cerr << "line " << line_num << ": ";
	throw xmlpp::exception (text);
    }
    
private:
    inline void confirm (bool condition, const Glib::ustring message) {
	if (!condition)
	    throw xmlpp::exception (message);
    }
    const Glib::ustring& get_attribute (const AttributeList& attributes, const Glib::ustring name, const Glib::ustring element) {
	AttributeList::const_iterator attribute = std::find_if(attributes.begin(), attributes.end(), AttributeHasName(name));
	if (attribute == attributes.end()) {
	    ostringstream msg;
	    msg << element << " element should have attribute \"" << name << '"';
	    throw xmlpp::exception (msg.str());
	}
	return attribute->value;
    }
    
    list<pair<string,TreeNode*> > elt_stack;
    map<string,TreeBranches*> refMap;	// for id / ref resolution of branches elements; contains pointers to all branches with an "id"
    size_t line_num;
};

int main(int argn, char** argv) {
    if (!(2 <= argn && argn <= 3)) {
	cerr << "Usage: "<<argv[0]<<" infile.xml [outfile.xml]" << endl;
	return 1;
    }
    const char *inFile = argv[1];
//     std::ostream& outStream = (argn == 3 ? ofstream(argv[2]) : cout);
    cout <<"Reading file "<< inFile<<endl;
    try {
	CMRefTreeParser parser;
	parser.set_substitute_entities (true);
	parser.parse_file (inFile);
    } catch (xmlpp::exception& e) {
	cerr << "error: " << e.what() << endl;
	return 1;
    }
    
    
    /*
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
    */
    return 0;
}
/*
uint recurseLeaves (Document!(char).Node node, ref double[uint] leaves, DecisionEnums inFlags = Decision::NONE, double inProb = 1.0) {
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
		    flags |= Decision::TEST_MICROSCOPY;
		} else if (end == "RDT") {
		    flags |= Decision::TEST_RDT;
		} else
		    Stdout.format ("unknown - {}", decis).newline;
	    } else if (prependedBy (decis, "result:", end)) {
		if (end == "true positive") {
		    flags |= Decision::RESULT_TRUE_POS;
		} else if (end == "false negative") {
		    flags |= Decision::RESULT_FALSE_NEG;
		} else if (end == "false positive") {
		    flags |= Decision::RESULT_FALSE_POS;
		} else if (end == "true negative") {
		    flags |= Decision::RESULT_TRUE_NEG;
		} else if (end == "positive") {
		    flags |= Decision::RESULT_FALSE;
		} else if (end == "negative") {
		    flags |= Decision::RESULT_TRUE;
		} else
		    Stdout.format ("unknown - {}", decis).newline;
	    } else if (prependedBy (decis, "drug:", end)) {
		if (end == "SP") {
		    flags |= Decision::DRUG_SP;
		} else if (end == "AL") {
		    flags |= Decision::DRUG_AL;
		} else if (end == "no antimalarial") {
		    flags |= Decision::DRUG_NO_AM;
		} else
		    Stdout.format ("unknown - {}", decis).newline;
	    } else if (prependedBy (decis, "quality:", end)) {
		if (end == "good") {
		    flags |= Decision::QUALITY_GOOD;
		} else if (end == "bad") {
		    flags |= Decision::QUALITY_BAD;
		} else
		    Stdout.format ("unknown - {}", decis).newline;
	    } else if (prependedBy (decis, "adherence:", end)) {
		if (end == "good") {
		    flags |= Decision::ADHERENCE_GOOD;
		} else if (end == "bad") {
		    flags |= Decision::ADHERENCE_BAD;
		} else
		    Stdout.format ("unknown - {}", decis).newline;
	    } else if (prependedBy (decis, "time:", end)) {
		int val = int.min;
		if (end.length == 1)
		    val = end[0] - '0';
		if (val >= 0 || val <= 15) {
		    flags |= (val << Decision::TSDELAY_LSHIFT);
		} else
		    Stdout.format ("unknown - {}", decis).newline;
	    } else if (prependedBy (decis, "treatment:", end)) {
		if (end == "no antimalarial") {
		    flags |= Decision::TREATMENT_NO_AM;
		} else if (end == "pre-referral only") {
		    flags |= Decision::TREATMENT_PREREF;
		} else {
		    flags |= Decision::DRUG_QN;
		    flags |= Decision::ADHERENCE_GOOD;
		    if (end == "parenteral") {
			flags |= Decision::TREATMENT_PARENTAL;
		    } else if (end == "delayed referral") {
			flags |= Decision::TREATMENT_DELREF | (1 << Decision::TSDELAY_LSHIFT);
		    } else if (end == "pre-referral and immediate referral") {
			flags |= Decision::TREATMENT_PREREF_IMMREF;
		    } else if (end == "pre-referral and delayed referral") {
			flags |= Decision::TREATMENT_PREREF_DELREF | (1 << Decision::TSDELAY_LSHIFT);
		    } else
			log.error ("unknown: {}", decis);
		}
	    } else if (prependedBy (decis, "management:", end)) {
		if (end == "good") {
		    flags |= Decision::MANAGEMENT_GOOD;
		} else if (end == "bad") {
		    flags |= Decision::MANAGEMENT_BAD;
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

bool prependedBy (string source, char[] match, out char[] end) {
    if (source.length() >= match.length() &&
	memcmp(source.c_str(), match.c_str(), match.length()) == 0) {
	end = Util.triml(source[match.length..$]);
	return true;
    }
    return false;
}

/** Take val & mask and right-shift by the number of binary zeros less-
 * significant than the first one in mask.
 * 
 * E.g. filterRShift!(0xF0) (0x40) == 0x4. * /
uint filterRShift(DecisionEnums mask, uint rShift = 0) (uint val) {
    static assert (rShift < uint.sizeof * 8);
    static if ((mask >> rShift) & 1) {
	return (val & mask) >> rShift;
    } else {
	return filterRShift!(mask, rShift+1) (val);
    }
}
*/