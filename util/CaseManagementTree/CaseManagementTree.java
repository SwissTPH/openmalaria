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

import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.Result;
import javax.xml.transform.Source;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

public class CM {

    DocumentBuilder _builder;
    
    String inFileName;
    String outFileName;

    public CM(String inN, String outN) throws ParserConfigurationException {
	inFileName = inN;
	outFileName = outN;
	_builder = DocumentBuilderFactory.newInstance()
                    .newDocumentBuilder();
    }

    public void processNodes(Element parent, Map<String, Double> outcomes,
            Double p) {
        Element childElement = (Element) parent.getChildNodes().item(1);
        if (childElement.getNodeName().equals("decision")) {

            String outcome = childElement.getTextContent();
            if (outcomes.containsKey(outcome)) {
                outcomes.put(outcome, outcomes.get(outcome) + p);
            } else {
                outcomes.put(outcome, p);
            }
            return;
        }

        NodeList nodeList = parent.getChildNodes();
        double pc = 0.0;
        for (int index = 0; index < nodeList.getLength(); index++) {
            if (nodeList.item(index).getNodeType() != Node.ELEMENT_NODE) {
                continue;
            }
            Element node = (Element) nodeList.item(index);
            Double pNode = Double.parseDouble(node.getAttribute("p"));
            processNodes(node, outcomes, p * pNode);
            pc = pc
                    + Double.parseDouble(((Element) nodeList.item(index))
                            .getAttribute("p"));
        }
        if (!(pc == 1.0)) {
            System.out.println("Parental node '"
                    + parent.getAttribute("decision")
                    + "' has at least one child with a wrong probability");
        }
        // System.out.println("prob. check of parent " +
        // parent.getAttribute("decision") + " = " + pc);
    }

    public double sumProbabilities(Element entryPoint) {
        double p = 0.0;
        NodeList decisionList = entryPoint.getChildNodes();
        for (int dcindex = 0; dcindex < decisionList.getLength(); dcindex++) {
            Element decision = (Element) decisionList.item(dcindex);
            p = p + Double.parseDouble(decision.getAttribute("p"));
        }
        return p;
    }

    public void checkUnitProbability() throws Exception {
	File outfile = new File (outFileName);
        Element caseManagements = _builder.parse(outfile).getDocumentElement();
        NodeList cmList = caseManagements
                .getElementsByTagName("caseManagement");
        for (int cmindex = 0; cmindex < cmList.getLength(); cmindex++) {
            NodeList cm = (NodeList) cmList.item(cmindex);
            for (int epindex = 0; epindex < cm.getLength(); epindex++) {
                Element ep = (Element) cm.item(epindex);
                if ((ep.getTagName() == "uc1") || (ep.getTagName() == "uc2")
                        || (ep.getTagName() == "sev")
                        || (ep.getTagName() == "nmfwithparasites")
                        || (ep.getTagName() == "nmfwithoutparasites")) {
                    double totalProb = sumProbabilities(ep);
                    if (!(totalProb == 1.0)) {
                        System.out
                                .println("Probabilities in Casemanagement number "
                                        + (cmindex + 1)
                                        + " entry point <"
                                        + ep.getTagName()
                                        + "> add up to "
                                        + totalProb + " instead of 1.0");
                    }
                }
            }
        }
    }

    public void translate() throws Exception {
	// source root:
	File file = new File (inFileName);
        Element ageDepDecisionTrees = _builder.parse(file).getDocumentElement();
	
	// target root:
        Document hsDoc = _builder.newDocument();
        Element newRoot = hsDoc.createElement("caseManagements");
        hsDoc.appendChild(newRoot);
        newRoot.setAttribute("name", "test");
	
        NodeList dtList = ageDepDecisionTrees
                .getElementsByTagName("decisionTrees");
        String minAY = "0.0";
        for (int dtindex = 0; dtindex < dtList.getLength(); dtindex++) {
            Element decisionTree = ((Element) dtList.item(dtindex));
            Element cm = hsDoc.createElement("caseManagement");
            newRoot.appendChild(cm);
            NodeList treeList = decisionTree.getElementsByTagName("tree");
            cm.setAttribute("minAgeYrs", minAY);
            cm
                    .setAttribute("maxAgeYrs", decisionTree
                            .getAttribute("maxAgeYrs"));
            minAY = decisionTree.getAttribute("maxAgeYrs").toString();
            for (int index = 0; index < treeList.getLength(); index++) {
                Map<String, Double> outcomeMap = new HashMap<String, Double>();
                Element tree = (Element) treeList.item(index);
                String entryPoint = tree.getAttribute("entryPoint");
                Element ep = hsDoc.createElement(entryPoint);
                cm.appendChild(ep);
                processNodes(tree, outcomeMap, 1.0);
                for (Iterator<String> outcomeiterator = outcomeMap.keySet()
                        .iterator(); outcomeiterator.hasNext();) {
                    String ocname = outcomeiterator.next();
                    Element oc = hsDoc.createElement("endPoint");
                    // oc.setAttribute("p",
                    // Double.toString(Math.round(outcomeMap
                    // .get(ocname) * 100.0) / 100.0));
                    oc.setAttribute("p", Double
                            .toString(outcomeMap.get(ocname)));
                    oc.setAttribute("decision", ocname);
                    ep.appendChild(oc);
                }
            }

            Node outcomes = ageDepDecisionTrees.getElementsByTagName(
                    "decisions").item(0);
            Element ocs = (Element) hsDoc.importNode(outcomes, true);
            cm.appendChild(ocs);
        }
        Source source = new DOMSource(hsDoc);
        // Prepare the output file
        File outFile = new File(outFileName);
        OutputStream os = new FileOutputStream(outFile);
        Result result = new StreamResult(os);
        // Write the DOM document to the file
        Transformer xformer = TransformerFactory.newInstance().newTransformer();
        xformer.transform(source, result);
    }

    public static int main(String[] args) {
	if (args.length != 2) {
	    System.err.println("Usage: "+args[0]+" IN.xml OUT.xml");
	    return 1;
	}
        try {
	    CM cm = new CM(args[0], args[1]);
	    cm.translate();
	    cm.checkUnitProbability();
	} catch (Exception e) {
            e.printStackTrace();
	    return -1;
        }
	
        System.out.println("\nDone");
	return 0;
    }
}
