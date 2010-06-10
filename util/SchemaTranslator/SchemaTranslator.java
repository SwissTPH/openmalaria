// This file is part of OpenMalaria.
// Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.StringReader;
import java.io.StringWriter;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Result;
import javax.xml.transform.Source;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.TransformerFactoryConfigurationError;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import javax.xml.validation.Schema;
import javax.xml.validation.SchemaFactory;
import javax.xml.validation.Validator;

import org.w3c.dom.Attr;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;

public class SchemaTranslator {

    DocumentBuilder _builder;
    Document scenarioDocument;
    Element scenarioElement;

    static final int CURRENT_VERSION = 19;

    private static int _required_version = CURRENT_VERSION;
    private static boolean doValidation = true;
    private static boolean doTranslation = true;
    private static boolean doODTTranslation = false;
    private static boolean doDBUpdate = false;

    private enum BugCorrectionBehaviour {
        none, correct, dontCorrect;
    }

    private static BugCorrectionBehaviour maxDensBug = BugCorrectionBehaviour.none;

    public static double HumanBloodIndex_NONNHS = 1;

    public static int INDEX_GAMBIAE_SS = 0;
    public static int INDEX_FUNESTUS = 1;
    public static int INDEX_ARABIENSIS = 2;

    public static String Name_GAMBIAE_SS = "gambiae_ss";
    public static String Name_FUNESTUS = "funestus";
    public static String Name_ARABIENSIS = "arabiensis";

    public static double[] HumanBloodIndexes = { 0.939, 0.98, 0.871 };
    public static double[] ProporitionsLaidEggsSameDay = { 0.313, 0.616, 0.313 };
    public static double[] PsSurvivalFeedingCycle = { 0.623, 0.611, 0.623 };
    public static double[] PAs = { 0.687, 0.384, 0.687 };
    public static double[] PA2s = { 0.0151, 0.00957, 0.320 };

    public static double td = 0.33;

    public static double Standard_RELATIVE_ENTO_AV = 1.0;
    public static double Standard_NHH_NUMBER = 1.0;

    public SchemaTranslator() {
        try {
            _builder = DocumentBuilderFactory.newInstance()
                    .newDocumentBuilder();
        } catch (ParserConfigurationException e) {
            e.printStackTrace();
        }
    }

    private void validate(Document scenarioDocument, String schemaFileName,
            String schemaDirectory) throws SAXException, IOException,
            Exception, TransformerFactoryConfigurationError {
        Document forValidation = (Document) scenarioDocument.cloneNode(true);
        Element scenarioElement = forValidation.getDocumentElement();
        scenarioElement.removeAttribute("xsi:noNamespaceSchemaLocation");

        // This is set by generateRun
        scenarioElement.setAttribute("assimMode", "0");
        // This is set by the work generator
        scenarioElement.setAttribute("wuID", "123");

        Element model = (Element) scenarioElement.getElementsByTagName("model")
                .item(0);
        Element t_parameters = (Element) model.getElementsByTagName(
                "parameters").item(0);

        if (t_parameters.getNodeValue() != null
                && t_parameters.getNodeValue().contains("@parameters@")) {
            t_parameters.getLastChild().setNodeValue("");
            // Add a dummy parameter
            Element parameters = forValidation.createElement("parameters");

            parameters.setAttribute("latentp", "0");
            parameters.setAttribute("delta", "0");
            parameters.setAttribute("interval", "0");
            parameters.setAttribute("iseed", "0");
            Element parameter = forValidation.createElement("parameter");
            parameter.setAttribute("value", "0");
            parameter.setAttribute("name", "0");
            parameter.setAttribute("number", "0");
            parameters.appendChild(parameter);
            scenarioElement.appendChild(parameters);
        }

        // Validate the updated document
        File xsdFile = new File(schemaDirectory + schemaFileName);
        if (xsdFile == null || !xsdFile.isFile()) {
            System.out.println("Unable to find " + schemaDirectory + schemaFileName
                    + " file; not validating.");
            return;
        }
        SchemaFactory factory = SchemaFactory
                .newInstance("http://www.w3.org/2001/XMLSchema");
        // NOTE: this may throw a java.lang.NullPointerException when run within
        // eclipse
        // run from a command-line instead (java SchemaTranslator).
        Schema schema = factory.newSchema(xsdFile);
        Validator validator = schema.newValidator();
        Source source = new DOMSource(forValidation);
        validator.validate(source);
    }
    
    // Helper function to strip old white-space.
    // Source: http://forums.sun.com/thread.jspa?threadID=5201482
    public static void visit(Node node, short nodeType, String name) {
	NodeList list = node.getChildNodes();
	for (int i = 0; i < list.getLength(); i++) {
	    // Get child node
	    Node childNode = list.item(i);
	    if (childNode.getNodeType() == nodeType && 
		(name == null || childNode.getNodeName().trim().equals(name) && childNode.getNodeValue().trim().equals(""))
	    ) {
		childNode.getParentNode().removeChild(childNode);
		// child was removed so list invalid; easiest is to start again:
		visit(node, nodeType, name);
		break;
	    }
	    else {
		visit(childNode, nodeType, name);
	    }
	}
    }
    private void translateFile(File documentFile, File outDir) throws Exception {
        scenarioDocument = _builder.parse(documentFile);
        String schemaFileName = translateDocument();
        if (schemaFileName == null) {
            System.err.println("Update of " + documentFile + " failed.");
            return;
        }
        if (doTranslation||doODTTranslation) {
            File outFile = new File(outDir, documentFile.getName());
            outFile.createNewFile();
            OutputStream os = new FileOutputStream(outFile);
            Result result = new StreamResult(os);
            // Write the DOM document to the file
            Transformer xformer = TransformerFactory.newInstance()
                    .newTransformer();
            xformer.setOutputProperty(OutputKeys.ENCODING, "UTF-8");
            xformer.setOutputProperty(OutputKeys.METHOD, "xml");
	    
	    // Reformat. First, remove all whitespace nodes, then insert new whitespace.
	    visit(scenarioElement, org.w3c.dom.Node.TEXT_NODE, "#text");
	    scenarioDocument.normalize();
	    // Then add new indentation/new-lines:
            xformer.setOutputProperty(OutputKeys.INDENT, "yes");
            xformer.setOutputProperty(
                    "{http://xml.apache.org/xslt}indent-amount", "2");
	    
            xformer.transform(new DOMSource(scenarioDocument), result);
        }
        if (doValidation)
            validate(scenarioDocument, schemaFileName, "../../test/");
    }

    private void updateDB() {
        Connection con = null;
        try {
            Class.forName("com.mysql.jdbc.Driver").newInstance();

            con = DriverManager.getConnection(
                    "jdbc:mysql://127.0.0.1:3306/DBNAME", "USER", "PASSWD");

            if (!con.isClosed())
                System.out.println("Successfully connected to "
                        + "MySQL server using TCP/IP...");
            Statement stmt = con.createStatement(ResultSet.TYPE_FORWARD_ONLY,
                    ResultSet.CONCUR_UPDATABLE);
            ResultSet uprs = stmt
                    .executeQuery("SELECT id,run_id,xml FROM scenarios");
            while (uprs.next()) {
                String xml = uprs.getString("xml");
                InputSource is = new InputSource();
                is.setCharacterStream(new StringReader(xml));
                _builder = DocumentBuilderFactory.newInstance()
                        .newDocumentBuilder();
                scenarioDocument = _builder.parse(is);
                String schemaFileName = translateDocument();
                Transformer transformer = TransformerFactory.newInstance()
                        .newTransformer();
                transformer.setOutputProperty(OutputKeys.INDENT, "yes");
                StreamResult result = new StreamResult(new StringWriter());
                DOMSource source = new DOMSource(scenarioDocument);
                transformer.transform(source, result);
                String translatedXML = result.getWriter().toString();
                uprs.updateString("xml", translatedXML);
                if (doValidation)
                    validate(scenarioDocument, schemaFileName, "./");
                uprs.updateRow();
            }
        } catch (SQLException e) {

            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            try {
                if (con != null)
                    con.close();
            } catch (SQLException e) {
            }
        }

    }

    private String translateDocument() throws NoSuchMethodException, Exception,
            IllegalAccessException, InvocationTargetException {
        scenarioElement = scenarioDocument.getDocumentElement();
        // System.out.println("Updating: " +
        // scenarioElement.getAttribute("name"));
        // 0 if no current version (getAttribute returns ""):
        int schemaVersion = Integer.parseInt("0"
                + scenarioElement.getAttribute("schemaVersion"));
        String schemaFileName = "scenario_" + schemaVersion + ".xsd";
        Class<? extends SchemaTranslator> cls = this.getClass();
        
        while (schemaVersion < _required_version) {
            ++schemaVersion;
            schemaFileName = "scenario_" + schemaVersion + ".xsd";
            scenarioElement.setAttribute("schemaVersion", Integer
                    .toString(schemaVersion));
            scenarioElement.setAttribute("xsi:noNamespaceSchemaLocation",
                    schemaFileName);

            String translateMeth = "translate" + (schemaVersion - 1) + "To"
                    + schemaVersion;
            Method method = cls.getMethod(translateMeth, new Class[] {});
            if (method == null)
                throw new Exception("Method " + translateMeth + " not found");
            if (!(Boolean) method.invoke(this, new Object[] {}))
                return null;
        }
        
        if(schemaVersion >= 18 && doODTTranslation)
        	if(!oDTTranslation())
        		return null;
        
        return schemaFileName;
    }
    
    

    // / Exactly what version 1 is has been forgotten; it's merged into 2.
    public Boolean translate0To1() {
        return true;
    }

    public Boolean translate1To2() {
        scenarioElement.setAttribute("xmlns:xsi",
                "http://www.w3.org/2001/XMLSchema-instance");
        // Done by translate() method
        // scenarioElement.setAttribute("xsi:noNamespaceSchemaLocation",
        // "scenario.xsd");
        if (!scenarioElement.hasAttribute("wuID"))
            scenarioElement.setAttribute("wuID", "0");
        if (!scenarioElement.hasAttribute("assimMode"))
            scenarioElement.setAttribute("assimMode", "0");
        Element elt = (Element) scenarioElement
                .getElementsByTagName("entoData").item(0);
        if (elt != null && elt.hasAttribute("firstDay")) {
            System.out.println("Warning: Removed firstDay attribute");
            elt.removeAttribute("firstDay");
        }
        elt = (Element) scenarioElement.getElementsByTagName("changeEIR").item(
                0);
        if (elt != null && elt.hasAttribute("firstDay")) {
            System.out.println("Warning: Removed firstDay attribute");
            elt.removeAttribute("firstDay");
        }
        elt = (Element) scenarioElement.getElementsByTagName("parameters")
                .item(0);
        if (elt != null && elt.hasAttribute("useIseed"))
            elt.removeAttribute("useIseed");

        NodeList sourcesElements = scenarioElement
                .getElementsByTagName("sources");
        for (int index = 0; index < sourcesElements.getLength(); index++) {
            Element sourcesElement = (Element) sourcesElements.item(index);
            Element parent = (Element) sourcesElement.getParentNode();
            parent.removeChild(sourcesElement);
            parent.setNodeValue("");
        }

        NodeList itemsElements = scenarioElement.getElementsByTagName("item");
        for (int index = 0; index < itemsElements.getLength(); index++) {
            Element itemsElement = (Element) itemsElements.item(index);
            itemsElement.setTextContent("");
        }
        Element paramElement = ((Element) scenarioElement.getElementsByTagName(
                "parameters").item(0));
        if (paramElement != null) {
            int nspore = Integer.parseInt(paramElement.getAttribute("nspore"));
            paramElement.setAttribute("eipDuration", Integer
                    .toString(nspore * 5));
            paramElement.removeAttribute("nspore");
        }
        NodeList mdaElements = scenarioElement.getElementsByTagName("MDA");
        for (int index = 0; index < mdaElements.getLength(); index++) {
            Element el = (Element) mdaElements.item(index);
            el.setAttribute("minAge", "0");
            el.setAttribute("maxAge", "99");
            el.setAttribute("coverage", "1");
        }
        NodeList allElements = scenarioElement.getElementsByTagName("*");
        for (int index = 0; index < allElements.getLength(); index++) {
            Element el = (Element) allElements.item(index);
            if (el.hasAttribute("best")) {
                String value = el.getAttribute("best");
                el.setAttribute("value", value);
                el.removeAttribute("best");
            }
        }
        return true;
    }

    public Boolean translate2To3() {
        return true;
    }

    /*
     * EntoData now has either a nonVector or a vector element; EIRDaily and
     * anopheles lists have moved to one of these. Some unwanted entomological
     * parameters have been removed, many have been added, and eipDuration has
     * been moved from param.
     */
    public Boolean translate3To4() {
        Element params = (Element) scenarioElement.getElementsByTagName(
                "parameters").item(0);
        Attr eip = null;
        if (params != null) {
            eip = params.getAttributeNode("eipDuration");
            params.removeAttributeNode(eip);
        }

        // Attribute added to nonVector, if used, below.

        Element elt = (Element) scenarioElement
                .getElementsByTagName("entoData").item(0);
        if (elt != null) {
            NodeList list = elt.getElementsByTagName("EIRDaily");
            if (list.getLength() != 0) {
                Element nonVector = scenarioDocument.createElement("nonVector");
                elt.appendChild(nonVector);
                // NOTE: results don't seem to be right if iterating forwards
                // instead of backwards
                for (int i = list.getLength() - 1; i >= 0; --i) {
                    // System.out.println (list.getLength());
                    Node eir = list.item(i);
                    // add the element at the new location, remove from old:
                    // unfortunately messes up the whitespace
                    nonVector.insertBefore(eir, nonVector.getFirstChild());
                }
                if (eip != null) {
                    nonVector.setAttributeNode(eip);
                } else {
                    nonVector.setAttribute("eipDuration", "10");
                }
            }
            list = elt.getElementsByTagName("anopheles");
            if (list.getLength() != 0) {
                Element vector = scenarioDocument.createElement("vector");
                elt.appendChild(vector);
                for (int i = list.getLength() - 1; i >= 0; --i) {
                    Node anoph = list.item(i);
                    ((Element) anoph).removeAttribute("useNv0Guess");
                    vector.insertBefore(anoph, vector.getFirstChild());
                }
            }

            elt.removeAttribute("inputType");
        }

        elt = (Element) scenarioElement.getElementsByTagName("interventions")
                .item(0);
        if (elt != null) {
            Element timed = (Element) elt.getElementsByTagName("timed").item(0);
            if (timed != null) {
                Element interv = (Element) timed.getElementsByTagName(
                        "intervention").item(0);
                if (interv != null) {
                    Element cEIR = (Element) interv.getElementsByTagName(
                            "changeEIR").item(0);
                    if (cEIR != null) {
                        cEIR.removeAttribute("inputType");
                        cEIR.removeAttribute("name"); // part of EntoData not
                        // NonVector
                        if (eip != null) {
                            cEIR.setAttributeNode((Attr) eip.cloneNode(true));
                        } else {
                            cEIR.setAttribute("eipDuration", "10");
                        }
                    }
                }
            }
        }
        return true;
    }

    // modelVersion flags 1<<2, 1<<4 or 1<<5 have changed
    public Boolean translate4To5() throws Exception {
        int ver = Integer
                .parseInt(scenarioElement.getAttribute("modelVersion"));
        if ((ver & 0x68) != 0) {// modelVersion with flags 1<<2, 1<<4 or 1<<5
            if ((ver & (1 << 5)) == 0) {
                ver = ver & (1 << 2);
            } else if ((ver & 0x68) == (1 << 5)) {
                ver = (ver ^ (1 << 5)) /* remove 1<<5 flag */
                        & (1 << 4); /* and add 1<<4 flag */
                System.err
                        .println("Warning: Scenario uses LOGNORMAL_MASS_ACTION_PLUS_PRE_IMM which has had a bug fixed!");
            } else {
                throw new Exception(
                        "Error: Scenario had a combination of InfectionIncidenceModel flags - this was invalid!");
            }
        }
        return true;
    }

    public Boolean translate5To6() throws Exception {
        int ver = Integer
                .parseInt(scenarioElement.getAttribute("modelVersion"));
        Element cMs = (Element) scenarioElement.getElementsByTagName(
                "caseManagements").item(0);
        // wuID is added by add_work.cpp
        scenarioElement.removeAttribute("wuID");
        if ((ver & 8192) != 0) { // ClinicalEventScheduler (new case
            // management)
            if (scenarioElement.getElementsByTagName("healthSystem")
                    .getLength() > 0)
                System.err
                        .println("Warning: healthSystem element present but not used");
        } else {
            if (cMs != null)
                System.err
                        .println("Warning: caseManagement element present but not used (updating anyway)");
        }
        if (cMs == null)
            return true; // element may not exist, in which case there's
        // nothing to
        // do
        NodeList cMList = cMs.getElementsByTagName("caseManagement");
        for (int i = 0; i < cMList.getLength(); ++i) {
            Element cM = (Element) cMList.item(i);
            cM.removeAttribute("minAgeYrs");
            Element nmfNP = (Element) cM.getElementsByTagName("nmf").item(0);
            scenarioDocument.renameNode(nmfNP, null, "nmfNP");
            Element nmfP = (Element) nmfNP.cloneNode(true);
            scenarioDocument.renameNode(nmfP, null, "nmfP");
            cM.insertBefore(nmfP, nmfNP);
        }
        return true;
    }

    // Version 7 added elements for ITN and IRS intervention descriptions.
    // Nothing old needs to be changed.
    public Boolean translate6To7() throws Exception {
        return true;
    }

    // Version 8 moved emergence rates and some other parameters into the XML
    // file. The relevant test scenarios have already been converted.
    public Boolean translate7To8() throws Exception {
        Element eD = (Element) scenarioElement.getElementsByTagName("entoData")
                .item(0);
        Element vect = (Element) eD.getElementsByTagName("vector").item(0);
        if (vect != null) {
            Element anoph = (Element) vect.getElementsByTagName("anopheles")
                    .item(0);
            Element mosq = (Element) anoph.getElementsByTagName("mosq").item(0);
            // This was required, so this if should always be true:
            if (mosq.getAttribute("emergenceRateFilename") != null) {
                System.err
                        .println("Warning: emergence rate data is now stored in the scenario document. Update by hand or run with \"openMalaria --enableERC\"");
                mosq.removeAttribute("emergenceRateFilename");
            }
        }
        return true;
    }

    // This changed some stuff to do with non-human hosts that wasn't used
    // before and added a VectorAvailability intervention.
    public Boolean translate8To9() throws Exception {
        return true;
    }

    // Version 10 introduced PKPD description parameters. No changes to
    // existing elements.
    public Boolean translate9To10() throws Exception {
        return true;
    }

    // Version 11 removes cached emerge rates from the schema
    public Boolean translate10To11() throws Exception {
        Element eD = (Element) scenarioElement.getElementsByTagName("entoData")
                .item(0);
        Element vect = (Element) eD.getElementsByTagName("vector").item(0);
        if (vect != null) {
            NodeList species = vect.getElementsByTagName("anopheles");
            for (int i = 0; i < species.getLength(); ++i) {
                Element anoph = (Element) species.item(i);
                Node er = anoph.getElementsByTagName("emergence").item(0);
                if (er != null)
                    anoph.removeChild(er);
                // These are from the parameter values based on Anopheles
                // gambiae in
                // Namawala, Tanzania, from the paper on comparing
                // interventions.
                anoph.setAttribute("propInfected", "0.078");
                anoph.setAttribute("propInfectious", "0.021");
            }
            System.err
                    .println("New attributes propInfected and propInfectious created with default values - please correct (for each anopheles section)!");
        }
        return true;
    }

    // Version 12 removed the simulationDuration attribute and changed the
    // event-scheduler data (no real scenarios yet so this is not auto-updated).
    public Boolean translate11To12() throws Exception {
        Element cms = (Element) scenarioElement.getElementsByTagName(
                "caseManagements").item(0);
        if (cms != null) {
            System.err
                    .println("Please replace the caseManagements element with an EventScheduler element (auto-update not implemented)");
            return false;
        }
        scenarioElement.removeAttribute("simulationDuration");
        return true;
    }

    // Version 13 replaced the modelVersion list with a ModelOptions section.
    // Similarly, the summaryOption attribute was replaced with a SurveyOptions
    // element.
    // Also, the GARKI_DENSITY_BIAS model option was introduced.
    // As such, old scenarios are definitely incompatible with the new code.
    public Boolean translate12To13() throws Exception {
        final int NUM_VERSIONS = 23;
        String[] num2String = new String[NUM_VERSIONS];
        num2String[1] = "PENALISATION_EPISODES";
        num2String[2] = "NEGATIVE_BINOMIAL_MASS_ACTION";
        num2String[3] = "ATTENUATION_ASEXUAL_DENSITY";
        num2String[4] = "LOGNORMAL_MASS_ACTION";
        num2String[5] = "NO_PRE_ERYTHROCYTIC";
        num2String[6] = "MAX_DENS_CORRECTION";
        num2String[7] = "INNATE_MAX_DENS";
        num2String[8] = "MAX_DENS_RESET";
        num2String[9] = "DUMMY_WITHIN_HOST_MODEL";
        num2String[10] = "PREDETERMINED_EPISODES";
        num2String[11] = "NON_MALARIA_FEVERS";
        num2String[12] = "INCLUDES_PK_PD";
        num2String[13] = "CLINICAL_EVENT_SCHEDULER";
        num2String[14] = "MUELLER_PRESENTATION_MODEL";
        num2String[15] = "TRANS_HET";
        num2String[16] = "COMORB_HET";
        num2String[17] = "TREAT_HET";
        num2String[18] = "COMORB_TRANS_HET";
        num2String[19] = "TRANS_TREAT_HET";
        num2String[20] = "COMORB_TREAT_HET";
        num2String[21] = "TRIPLE_HET";
        num2String[22] = "EMPIRICAL_WITHIN_HOST_MODEL";

        Element modelOptions = scenarioDocument.createElement("ModelOptions");
        int verFlags = Integer.parseInt(scenarioElement
                .getAttribute("modelVersion"));

        for (int i = 1; i < NUM_VERSIONS; ++i) {
            if ((verFlags & (1 << i)) != 0) {
                Element opt = scenarioDocument.createElement("option");
                opt.setAttribute("name", num2String[i]);
                opt.setAttribute("value", "true");
                modelOptions.appendChild(opt);
            }
        }
        if ((verFlags & (1 << 6)) == 0 && (verFlags & (1 << 9)) == 0
                && (verFlags & (1 << 22)) == 0) {
            if (maxDensBug == BugCorrectionBehaviour.correct) {
                // option is enabled by default so don't need to add it
            } else if (maxDensBug == BugCorrectionBehaviour.dontCorrect) {
                Element opt = scenarioDocument.createElement("option");
                opt.setAttribute("name", num2String[6]);
                opt.setAttribute("value", "false");
                modelOptions.appendChild(opt);
            } else {
                System.err
                        .println("scenario doesn't include MAX_DENS_CORRECTION: please specify --maxDensCorrection BOOL");
                return false;
            }
        }

        scenarioElement.insertBefore(modelOptions, scenarioElement
                .getFirstChild());
        scenarioElement.removeAttribute("modelVersion");

        final int NUM_SURVEYS = 31;
        String[] surveyCode2String = new String[NUM_SURVEYS];
        surveyCode2String[0] = "nHost";
        surveyCode2String[1] = "nInfect";
        surveyCode2String[2] = "nExpectd";
        surveyCode2String[3] = "nPatent";
        surveyCode2String[4] = "sumLogPyrogenThres";
        surveyCode2String[5] = "sumlogDens";
        surveyCode2String[6] = "totalInfs";
        surveyCode2String[7] = "nTransmit";
        surveyCode2String[8] = "totalPatentInf";
        surveyCode2String[9] = "contrib";
        surveyCode2String[10] = "sumPyrogenThresh";
        surveyCode2String[11] = "nTreatments1";
        surveyCode2String[12] = "nTreatments2";
        surveyCode2String[13] = "nTreatments3";
        surveyCode2String[14] = "nUncomp";
        surveyCode2String[15] = "nSevere";
        surveyCode2String[16] = "nSeq";
        surveyCode2String[17] = "nHospitalDeaths";
        surveyCode2String[18] = "nIndDeaths";
        surveyCode2String[19] = "nDirDeaths";
        surveyCode2String[20] = "nEPIVaccinations";
        surveyCode2String[21] = "imr_summary";
        surveyCode2String[22] = "nMassVaccinations";
        surveyCode2String[23] = "nHospitalRecovs";
        surveyCode2String[24] = "nHospitalSeqs";
        surveyCode2String[25] = "nIPTDoses";
        surveyCode2String[26] = "annAvgK";
        surveyCode2String[27] = "nNMFever";
        surveyCode2String[28] = "innoculationsPerDayOfYear";
        surveyCode2String[29] = "kappaPerDayOfYear";
        surveyCode2String[30] = "innoculationsPerAgeGroup";

        Element surveyOptions = scenarioDocument.createElement("SurveyOptions");
        Element monitoring = (Element) scenarioElement.getElementsByTagName(
                "monitoring").item(0);
        Element surveys = (Element) monitoring.getElementsByTagName("surveys")
                .item(0);
        int surveyFlags = Integer.parseInt(surveys
                .getAttribute("summaryOption"));

        for (int i = 0; i < NUM_SURVEYS; ++i) {
            if ((surveyFlags & (1 << i)) != 0) {
                Element opt = scenarioDocument.createElement("option");
                opt.setAttribute("name", surveyCode2String[i]);
                opt.setAttribute("value", "true");
                surveyOptions.appendChild(opt);
            }
        }

        monitoring.insertBefore(surveyOptions, surveys);
        surveys.removeAttribute("summaryOption");

        int analysisNum = Integer.parseInt(scenarioElement
                .getAttribute("analysisNo"));
        // These analysis numbers _were_ reserved for Garki scenarios.
        if ((analysisNum >= 22) && (analysisNum <= 30)) {
            Element opt = scenarioDocument.createElement("option");
            opt.setAttribute("name", "GARKI_DENSITY_BIAS");
            opt.setAttribute("value", "true");
            modelOptions.appendChild(opt);
        }

        return true;
    }

    // Version 14 changed the drugDescription element. This was as-yet unused
    // and there's no direct
    // translation from the old version.
    public Boolean translate13To14() throws Exception {
        Element cms = (Element) scenarioElement.getElementsByTagName(
                "drugDescription").item(0);
        if (cms != null) {
            System.err
                    .println("Warning: drugDescription element has changed; please rewrite manually.");
        }
        return true;
    }

    // Version 15 allowed MDA interventions to include drug information (no
    // changes to existing scenarios)
    public Boolean translate14To15() throws Exception {
        return true;
    }

    public Boolean translate15To16() throws Exception {

        Element model = scenarioDocument.createElement("model");
        Element clinical = scenarioDocument.createElement("clinical");
        Element modelOptions = (Element) scenarioElement.getElementsByTagName(
                "ModelOptions").item(0);
        Element parameters = (Element) scenarioElement.getElementsByTagName(
                "parameters").item(0);

        model.appendChild(modelOptions);
        model.appendChild(clinical);
        model.appendChild(parameters);

        scenarioElement.appendChild(model);

        Element healthSystemOld = (Element) scenarioElement
                .getElementsByTagName("healthSystem").item(0);
        Element eventScheduler = (Element) scenarioElement
                .getElementsByTagName("EventScheduler").item(0);
        Attr healthSystemMemory;

        Element healthSystemNew = scenarioDocument
                .createElement("healthSystem");

        if (healthSystemOld == null) {
            healthSystemMemory = eventScheduler
                    .getAttributeNode("healthSystemMemory");
            eventScheduler.removeAttribute("healthSystemMemory");
            // scenarioDocument.renameNode(eventScheduler, null,
            // "HSEventScheduler");
            Element CFR = scenarioDocument.createElement("CFR");
            Element group = scenarioDocument.createElement("group");
            Attr cfr = scenarioDocument.createAttribute("cfr");
            Attr lowerbound = scenarioDocument.createAttribute("lowerbound");

            cfr.setNodeValue("0");
            lowerbound.setNodeValue("0");

            group.setAttributeNode(cfr);
            group.setAttributeNode(lowerbound);

            CFR.appendChild(group);

            healthSystemNew.appendChild(eventScheduler);
            healthSystemNew.appendChild(CFR);
        } else {
            healthSystemMemory = healthSystemOld
                    .getAttributeNode("healthSystemMemory");
            healthSystemOld.removeAttribute("healthSystemMemory");
            // healthSystemOld.removeAttribute("name");
            scenarioDocument.renameNode(healthSystemOld, null,
                    "ImmediateOutcomes");
            Element CFR = (Element) healthSystemOld.getElementsByTagName("CFR")
                    .item(0);

            healthSystemNew.appendChild(healthSystemOld);
            healthSystemNew.appendChild(CFR);
        }

        Element Intervention = (Element) scenarioElement.getElementsByTagName(
                "intervention").item(0);
        if (Intervention != null) {
            Element changeHS = (Element) Intervention.getElementsByTagName(
                    "changeHS").item(0);

            if (changeHS != null) {
                changeHS.removeAttribute("healthSystemMemory");

                scenarioDocument
                        .renameNode(changeHS, null, "ImmediateOutcomes");

                Element changeHSNew = scenarioDocument
                        .createElement("changeHS");
                Intervention.appendChild(changeHSNew);
                changeHSNew.appendChild(changeHS);

                Element HSCFR = (Element) changeHS.getElementsByTagName("CFR")
                        .item(0);
                changeHSNew.appendChild(HSCFR);
            }

        }

        scenarioElement.insertBefore(healthSystemNew, scenarioElement
                .getElementsByTagName("entoData").item(0));
        clinical.setAttributeNode(healthSystemMemory);

        return true;
    }

    public boolean translate16To17() throws Exception {

        Element vector = (Element) scenarioElement.getElementsByTagName(
                "vector").item(0);

        if (vector != null) {
            NodeList anopheles = vector.getElementsByTagName("anopheles");

            if (((Element) anopheles.item(0)).getElementsByTagName(
                    "nonHumanHosts").getLength() > 0) {
                Element nhh = (Element) ((Element) anopheles.item(0))
                        .getElementsByTagName("nonHumanHosts").item(0);

                Attr nonHumanHostsnumber = scenarioDocument
                        .createAttribute("number");
                nonHumanHostsnumber.setNodeValue(Double
                        .toString(SchemaTranslator.Standard_NHH_NUMBER));
                Attr name = scenarioDocument.createAttribute("name");
                name.setNodeValue(nhh.getAttribute("name"));

                Element nonHumanHostsnumbers = scenarioDocument
                        .createElement("nonHumanHosts");
                nonHumanHostsnumbers.setAttributeNode(name);
                nonHumanHostsnumbers.setAttributeNode(nonHumanHostsnumber);

                vector.appendChild(nonHumanHostsnumbers);
            }

            for (int i = 0; i < anopheles.getLength(); i++) {
                Element anophelesType = (Element) anopheles.item(i);
                String typeName = anophelesType.getAttribute("mosquito");
                NodeList nonHumanHosts = anophelesType
                        .getElementsByTagName("nonHumanHosts");
                Element mosq = (Element) anophelesType.getElementsByTagName(
                        "mosq").item(0);

                if (typeName.equals(SchemaTranslator.Name_GAMBIAE_SS))
                    setMosqsNewAttributes(SchemaTranslator.INDEX_GAMBIAE_SS,
                            mosq, nonHumanHosts);

                else if (typeName.equals(SchemaTranslator.Name_FUNESTUS))
                    setMosqsNewAttributes(SchemaTranslator.INDEX_FUNESTUS,
                            mosq, nonHumanHosts);

                else if (typeName.equals(SchemaTranslator.Name_ARABIENSIS))
                    setMosqsNewAttributes(SchemaTranslator.INDEX_ARABIENSIS,
                            mosq, nonHumanHosts);

                else {
                    System.err
                            .println("There are no standards values for this kind of mosquito. Please edit those values per hand. ");
                    System.err.println("This scenario will not be updated.");
                    return false;
                }
            }
        }
        return true;
    }

    public boolean translate17To18() throws Exception {

        Attr popSize = scenarioElement.getAttributeNode("popSize");
        Attr maxAgeYrs = scenarioElement.getAttributeNode("maximumAgeYrs");
        Attr mode = scenarioElement.getAttributeNode("mode");

        scenarioElement.removeAttribute("popSize");
        scenarioElement.removeAttribute("maximumAgeYrs");
        scenarioElement.removeAttribute("mode");

        Element demography = (Element) scenarioElement.getElementsByTagName(
                "demography").item(0);
        Element entoData = (Element) scenarioElement.getElementsByTagName(
                "entoData").item(0);

        demography.setAttributeNode(popSize);
        demography.setAttributeNode(maxAgeYrs);

        entoData.setAttributeNode(mode);

        return true;
    }
    
    // Version 19: entoData's "mode" can no longer have value 3.
    // Removed unused "delta" from parameters.
    // Moved two event-scheduler outcome attributes into parameters element -- updated by hand.
    public boolean translate18To19() throws Exception {
	Element entoData = (Element) scenarioElement.getElementsByTagName(
	    "entoData").item(0);
	Attr mode = entoData.getAttributeNode("mode");
	
	if (Integer.parseInt(mode.getValue()) == 3){
	    // unless an intervention at time 0 specifies EIR values, the scenario was buggy
	    boolean hasTransientEIRAt0 = false;
	    
	    Element elt = (Element) scenarioElement.getElementsByTagName("interventions").item(0);
	    if (elt != null) {
		Element timed = (Element) elt.getElementsByTagName("timed").item(0);
		if (timed != null) {
		    NodeList intervs = timed.getElementsByTagName("intervention");
		    for (int i = 0; i < intervs.getLength(); i++) {
			Element interv = (Element) intervs.item(i);
			if (Integer.parseInt(interv.getAttribute("time")) != 0)
			    continue;	// only interested in interv. at time 0
			
			Element cEIR = (Element) interv.getElementsByTagName("changeEIR").item(0);
			if (cEIR != null)	// yes, have applicable transient EIR data
			    hasTransientEIRAt0 = true;
		    }
		}
	    }
	    if (!hasTransientEIRAt0) {
		System.err.println("Error: entoData has mode=\"3\", but no changeEIR intervention found at time 0.");
		return false;
	    }
	    
	    // Assuming correct changeEIR intervention was found, we can just switch mode to 4.
	    mode.setValue ("4");
	}
	
        Element model = (Element) scenarioElement.getElementsByTagName("model")
                .item(0);
        if(model!=null)
        {
            Element t_parameters = (Element) model.getElementsByTagName(
                "parameters").item(0);
            if(t_parameters!=null)
                t_parameters.removeAttribute("delta");
        }
        
	
	return true;
    }

    /**
     * This function is used to translate the 5-day timestep fitting
     * scenarii to 1-day timestep fitting scenarii. Since we're using a fairly
     * simple case management description (no interventions, no treatment),
     * then it's not too difficult to translate those scenarii. This translation
     * is therefore not intended for more complicated 5-day timestep scenarii.
     * @return true if the translation was a success.
     *
     */
    private Boolean oDTTranslation() {

    	Element surveys = (Element)scenarioElement.getElementsByTagName("surveys").item(0);
        NodeList surveystimes = surveys.getElementsByTagName("surveyTime");

        for(int i=0;i<surveystimes.getLength();i++)
        {
            Element surveytime = (Element)surveystimes.item(i);
            int surveyTime = Integer.parseInt(surveytime.getTextContent());
            surveytime.setTextContent(String.valueOf(((surveyTime -1)*5)+1));
        }

        Element modelElement = (Element)scenarioElement.getElementsByTagName("model").item(0);
        Element modelOptions = (Element)modelElement.getElementsByTagName("ModelOptions").item(0);

        Element molineauxOption = scenarioDocument.createElement("option");
        Attr name = scenarioDocument.createAttribute("name");
        name.setNodeValue("MOLINEAUX_WITHIN_HOST_MODEL");
        Attr value = scenarioDocument.createAttribute("value");
        value.setNodeValue("true");
        molineauxOption.setAttributeNode(name);
        molineauxOption.setAttributeNode(value);
        modelOptions.appendChild(molineauxOption);

        Element pkpdOption = scenarioDocument.createElement("option");
        Attr namePkPd = scenarioDocument.createAttribute("name");
        namePkPd.setNodeValue("INCLUDES_PK_PD");
        Attr valuePkPd = scenarioDocument.createAttribute("value");
        valuePkPd.setNodeValue("true");
        pkpdOption.setAttributeNode(namePkPd);
        pkpdOption.setAttributeNode(valuePkPd);
        modelOptions.appendChild(pkpdOption);

        Element esOption = scenarioDocument.createElement("option");
        Attr nameES = scenarioDocument.createAttribute("name");
        nameES.setNodeValue("CLINICAL_EVENT_SCHEDULER");
        Attr valueES = scenarioDocument.createAttribute("value");
        valueES.setNodeValue("true");
        esOption.setAttributeNode(nameES);
        esOption.setAttributeNode(valueES);
        modelOptions.appendChild(esOption);

        Element clinical = (Element)scenarioElement.getElementsByTagName("clinical").item(0);
        Attr healthSystemMemory = clinical.getAttributeNode("healthSystemMemory");
        healthSystemMemory.setValue(String.valueOf(28));

        NodeList interventionList = (NodeList)scenarioElement.getElementsByTagName("intervention");

        for(int i=0;i<interventionList.getLength();i++)
        {
            Element intervention = (Element)interventionList.item(0);
            Attr time = intervention.getAttributeNode("time");

            if(time!=null)
                time.setNodeValue(String.valueOf(((Integer.parseInt(time.getNodeValue())-1)*5)+1));
        }

        NodeList changeHSList = (NodeList)scenarioElement.getElementsByTagName("changeHS");

        for(int i=0;i<changeHSList.getLength();i++)
        {
            Element changeHS = (Element) changeHSList.item(i);
            Element immediateOutcomes = (Element)changeHS.getElementsByTagName("ImmediateOutcomes").item(0);
            String valueString = immediateOutcomes.getAttribute("name");

            if(valueString.equals("Do Monitoring HS")||valueString.equals("Np Monitoring HS"))
            {
                Element eventScheduler = (Element)scenarioDocument.createElement("EventScheduler");
                changeHS.removeChild(immediateOutcomes);

                Element uncomplicated = (Element)scenarioDocument.createElement("uncomplicated");
                Element decisions = (Element)scenarioDocument.createElement("decisions");

                Element decisionTreat = (Element)scenarioDocument.createElement("decision");

                Attr nameTreat = scenarioDocument.createAttribute("name");
                nameTreat.setNodeValue("treatment");
                Attr dependsTreat = scenarioDocument.createAttribute("depends");
                dependsTreat.setNodeValue("");
                Attr valuesTreat = scenarioDocument.createAttribute("values");
                valuesTreat.setNodeValue("effective_treat,none");

                decisionTreat.setAttributeNode(nameTreat);
                decisionTreat.setAttributeNode(dependsTreat);
                decisionTreat.setAttributeNode(valuesTreat);
                decisionTreat.setTextContent("effective_treat");
                decisions.appendChild(decisionTreat);

                Element decisionTest = (Element)scenarioDocument.createElement("decision");

                Attr nameTest = scenarioDocument.createAttribute("name");
                nameTest.setNodeValue("test");
                Attr dependsTest = scenarioDocument.createAttribute("depends");
                dependsTest.setNodeValue("");
                Attr valuesTest = scenarioDocument.createAttribute("values");
                valuesTest.setNodeValue("none,microscopy,RDT");

                decisionTest.setAttributeNode(nameTest);
                decisionTest.setAttributeNode(dependsTest);
                decisionTest.setAttributeNode(valuesTest);
                decisionTest.setTextContent("none");
                decisions.appendChild(decisionTest);

                uncomplicated.appendChild(decisions);

                Element treatments = (Element)scenarioDocument.createElement("treatments");
                Element treatment = (Element)scenarioDocument.createElement("treatment");
                Attr nameTreatEl = scenarioDocument.createAttribute("name");
                nameTreatEl.setNodeValue("effective_treat");
                treatment.setAttributeNode(nameTreatEl);

                Element schedule = (Element)scenarioDocument.createElement("schedule");
                Element medicate = (Element)scenarioDocument.createElement("medicate");

                Attr drug = (Attr)scenarioDocument.createAttribute("drug");
                drug.setNodeValue("effective");
                Attr mg = (Attr)scenarioDocument.createAttribute("mg");
                mg.setNodeValue("1");
                Attr hour = (Attr)scenarioDocument.createAttribute("hour");
                hour.setNodeValue("0");

                medicate.setAttributeNode(drug);
                medicate.setAttributeNode(mg);
                medicate.setAttributeNode(hour);

                schedule.appendChild(medicate);
                treatment.appendChild(schedule);
                treatments.appendChild(treatment);

                Element treatmentNone = (Element)scenarioDocument.createElement("treatment");
                Attr nameTreatNone = (Attr)scenarioDocument.createAttribute("name");
                nameTreatNone.setNodeValue("none");
                treatmentNone.setAttributeNode(nameTreatNone);

                Element scheduleNone = (Element)scenarioDocument.createElement("schedule");
                treatmentNone.appendChild(scheduleNone);
                treatments.appendChild(treatmentNone);
                uncomplicated.appendChild(treatments);

                eventScheduler.appendChild(uncomplicated);

                Element complicated = (Element)scenarioDocument.createElement("complicated");
                Element decisionsComp = (Element)scenarioDocument.createElement("decisions");

                Element decisionCompTreat = (Element)scenarioDocument.createElement("decision");

                Attr nameCompTreat = scenarioDocument.createAttribute("name");
                nameCompTreat.setNodeValue("treatment");
                Attr dependsCompTreat = scenarioDocument.createAttribute("depends");
                dependsCompTreat.setNodeValue("");
                Attr valuesCompTreat = scenarioDocument.createAttribute("values");
                valuesCompTreat.setNodeValue("effective_treat,none");

                decisionCompTreat.setAttributeNode(nameCompTreat);
                decisionCompTreat.setAttributeNode(dependsCompTreat);
                decisionCompTreat.setAttributeNode(valuesCompTreat);
                decisionCompTreat.setTextContent("effective_treat");
                decisionsComp.appendChild(decisionCompTreat);

                Element decisionCompHosp = (Element)scenarioDocument.createElement("decision");

                Attr nameCompHosp = scenarioDocument.createAttribute("name");
                nameCompHosp.setNodeValue("hospitalisation");
                Attr dependsCompHosp = scenarioDocument.createAttribute("depends");
                dependsCompHosp.setNodeValue("");
                Attr valuesCompHosp = scenarioDocument.createAttribute("values");
                valuesCompHosp.setNodeValue("none,delayed,immediate");

                decisionCompHosp.setAttributeNode(nameCompHosp);
                decisionCompHosp.setAttributeNode(dependsCompHosp);
                decisionCompHosp.setAttributeNode(valuesCompHosp);
                decisionCompHosp.setTextContent("immediate");
                decisionsComp.appendChild(decisionCompHosp);

                Element decisionCompTest = (Element)scenarioDocument.createElement("decision");

                Attr nameCompTest = scenarioDocument.createAttribute("name");
                nameCompTest.setNodeValue("test");
                Attr dependsCompTest = scenarioDocument.createAttribute("depends");
                dependsCompTest.setNodeValue("");
                Attr valuesCompTest = scenarioDocument.createAttribute("values");
                valuesCompTest.setNodeValue("none,microscopy,RDT");

                decisionCompTest.setAttributeNode(nameCompTest);
                decisionCompTest.setAttributeNode(dependsCompTest);
                decisionCompTest.setAttributeNode(valuesCompTest);
                decisionCompTest.setTextContent("none");
                decisionsComp.appendChild(decisionCompTest);

                complicated.appendChild(decisionsComp);

                Element treatmentsComp = (Element)scenarioDocument.createElement("treatments");
                Element treatmentComp = (Element)scenarioDocument.createElement("treatment");
                Attr nameTreatComp = scenarioDocument.createAttribute("name");
                nameTreatComp.setNodeValue("effective_treat");
                treatmentComp.setAttributeNode(nameTreatComp);

                Element scheduleComp = (Element)scenarioDocument.createElement("schedule");
                Element medicateComp = (Element)scenarioDocument.createElement("medicate");

                Attr drugComp = (Attr)scenarioDocument.createAttribute("drug");
                drugComp.setNodeValue("effective");
                Attr mgComp = (Attr)scenarioDocument.createAttribute("mg");
                mgComp.setNodeValue("1");
                Attr hourComp = (Attr)scenarioDocument.createAttribute("hour");
                hourComp.setNodeValue("0");

                medicateComp.setAttributeNode(drugComp);
                medicateComp.setAttributeNode(mgComp);
                medicateComp.setAttributeNode(hourComp);

                scheduleComp.appendChild(medicateComp);
                treatmentComp.appendChild(scheduleComp);
                treatmentsComp.appendChild(treatmentComp);

                Element treatmentCompNone = (Element)scenarioDocument.createElement("treatment");
                Attr nameTreatCompNone = (Attr)scenarioDocument.createAttribute("name");
                nameTreatCompNone.setNodeValue("none");
                treatmentCompNone.setAttributeNode(nameTreatCompNone);

                Element scheduleCompNone = (Element)scenarioDocument.createElement("schedule");
                treatmentCompNone.appendChild(scheduleCompNone);
                treatmentsComp.appendChild(treatmentCompNone);
                complicated.appendChild(treatmentsComp);

                eventScheduler.appendChild(complicated);

                Element clinicalOutcomes = (Element)scenarioDocument.createElement("ClinicalOutcomes");
                Element maxUCSeekingMemory = (Element)scenarioDocument.createElement("maxUCSeekingMemory");
                maxUCSeekingMemory.setTextContent("3");
                Element uncomplicatedCaseDuration = (Element)scenarioDocument.createElement("uncomplicatedCaseDuration");
                uncomplicatedCaseDuration.setTextContent("3");
                Element complicatedCaseDuration = (Element)scenarioDocument.createElement("complicatedCaseDuration");
                complicatedCaseDuration.setTextContent("5");
                Element complicatedRiskDuration = (Element)scenarioDocument.createElement("complicatedRiskDuration");
                complicatedRiskDuration.setTextContent("5");
                Element pImmediateUC = (Element)scenarioDocument.createElement("pImmediateUC");
                pImmediateUC.setTextContent("1");
                //Element propDeathsFirstDay = (Element)scenarioDocument.createElement("propDeathsFirstDay");
                //propDeathsFirstDay.setTextContent("0.4");

                //this communityOddsMultiplier will be removed (for schema >= 19)
                //Element communityOddsMultiplier = (Element)scenarioDocument.createElement("communityOddsMultiplier");
                //communityOddsMultiplier.setTextContent("1.5");

                clinicalOutcomes.appendChild(maxUCSeekingMemory);
                clinicalOutcomes.appendChild(uncomplicatedCaseDuration);
                clinicalOutcomes.appendChild(complicatedCaseDuration);
                clinicalOutcomes.appendChild(complicatedRiskDuration);
                clinicalOutcomes.appendChild(pImmediateUC);
                

                //clinicalOutcomes.appendChild(communityOddsMultiplier);

                eventScheduler.appendChild(clinicalOutcomes);


                Element CFR = (Element)changeHS.getElementsByTagName("CFR").item(0);
                if(CFR!=null)
                    changeHS.insertBefore(eventScheduler, CFR);
                else
                    changeHS.appendChild(eventScheduler);
            }
            else return false;
        }



        Element healthSystem = (Element)scenarioDocument.getElementsByTagName("healthSystem").item(0);

        Element immediateOutcomes = (Element)healthSystem.getElementsByTagName("ImmediateOutcomes").item(0);
        String valueString = immediateOutcomes.getAttribute("name");

        System.out.println(valueString);

        if(valueString.equals("no Treatment")||valueString.equals("Mortality Fitting")||valueString.equals("no Treatment no Mortality"))
        {
            Element eventScheduler = (Element)scenarioDocument.createElement("EventScheduler");
            healthSystem.removeChild(immediateOutcomes);
            
            Element uncomplicated = (Element)scenarioDocument.createElement("uncomplicated");
            Element decisions = (Element)scenarioDocument.createElement("decisions");

            Element decisionTreat = (Element)scenarioDocument.createElement("decision");

            Attr nameTreat = scenarioDocument.createAttribute("name");
            nameTreat.setNodeValue("treatment");
            Attr dependsTreat = scenarioDocument.createAttribute("depends");
            dependsTreat.setNodeValue("");
            Attr valuesTreat = scenarioDocument.createAttribute("values");
            valuesTreat.setNodeValue("effective_treat,none");

            decisionTreat.setAttributeNode(nameTreat);
            decisionTreat.setAttributeNode(dependsTreat);
            decisionTreat.setAttributeNode(valuesTreat);
            decisionTreat.setTextContent("none");
            decisions.appendChild(decisionTreat);

            Element decisionTest = (Element)scenarioDocument.createElement("decision");

            Attr nameTest = scenarioDocument.createAttribute("name");
            nameTest.setNodeValue("test");
            Attr dependsTest = scenarioDocument.createAttribute("depends");
            dependsTest.setNodeValue("");
            Attr valuesTest = scenarioDocument.createAttribute("values");
            valuesTest.setNodeValue("none,microscopy,RDT");

            decisionTest.setAttributeNode(nameTest);
            decisionTest.setAttributeNode(dependsTest);
            decisionTest.setAttributeNode(valuesTest);
            decisionTest.setTextContent("none");
            decisions.appendChild(decisionTest);

            uncomplicated.appendChild(decisions);

            Element treatments = (Element)scenarioDocument.createElement("treatments");
            Element treatment = (Element)scenarioDocument.createElement("treatment");
            Attr nameTreatEl = scenarioDocument.createAttribute("name");
            nameTreatEl.setNodeValue("effective_treat");
            treatment.setAttributeNode(nameTreatEl);

            Element schedule = (Element)scenarioDocument.createElement("schedule");
            Element medicate = (Element)scenarioDocument.createElement("medicate");

            Attr drug = (Attr)scenarioDocument.createAttribute("drug");
            drug.setNodeValue("effective");
            Attr mg = (Attr)scenarioDocument.createAttribute("mg");
            mg.setNodeValue("1");
            Attr hour = (Attr)scenarioDocument.createAttribute("hour");
            hour.setNodeValue("0");

            medicate.setAttributeNode(drug);
            medicate.setAttributeNode(mg);
            medicate.setAttributeNode(hour);

            schedule.appendChild(medicate);
            treatment.appendChild(schedule);
            treatments.appendChild(treatment);

            Element treatmentNone = (Element)scenarioDocument.createElement("treatment");
            Attr nameTreatNone = (Attr)scenarioDocument.createAttribute("name");
            nameTreatNone.setNodeValue("none");
            treatmentNone.setAttributeNode(nameTreatNone);

            Element scheduleNone = (Element)scenarioDocument.createElement("schedule");
            treatmentNone.appendChild(scheduleNone);
            treatments.appendChild(treatmentNone);
            uncomplicated.appendChild(treatments);
            
            eventScheduler.appendChild(uncomplicated);

            Element complicated = (Element)scenarioDocument.createElement("complicated");
            Element decisionsComp = (Element)scenarioDocument.createElement("decisions");

            Element decisionCompTreat = (Element)scenarioDocument.createElement("decision");

            Attr nameCompTreat = scenarioDocument.createAttribute("name");
            nameCompTreat.setNodeValue("treatment");
            Attr dependsCompTreat = scenarioDocument.createAttribute("depends");
            dependsCompTreat.setNodeValue("");
            Attr valuesCompTreat = scenarioDocument.createAttribute("values");
            valuesCompTreat.setNodeValue("effective_treat,none");

            decisionCompTreat.setAttributeNode(nameCompTreat);
            decisionCompTreat.setAttributeNode(dependsCompTreat);
            decisionCompTreat.setAttributeNode(valuesCompTreat);
            decisionCompTreat.setTextContent("none");
            decisionsComp.appendChild(decisionCompTreat);

            Element decisionCompHosp = (Element)scenarioDocument.createElement("decision");

            Attr nameCompHosp = scenarioDocument.createAttribute("name");
            nameCompHosp.setNodeValue("hospitalisation");
            Attr dependsCompHosp = scenarioDocument.createAttribute("depends");
            dependsCompHosp.setNodeValue("");
            Attr valuesCompHosp = scenarioDocument.createAttribute("values");
            valuesCompHosp.setNodeValue("none,delayed,immediate");

            decisionCompHosp.setAttributeNode(nameCompHosp);
            decisionCompHosp.setAttributeNode(dependsCompHosp);
            decisionCompHosp.setAttributeNode(valuesCompHosp);
            decisionCompHosp.setTextContent("none");
            decisionsComp.appendChild(decisionCompHosp);

            Element decisionCompTest = (Element)scenarioDocument.createElement("decision");

            Attr nameCompTest = scenarioDocument.createAttribute("name");
            nameCompTest.setNodeValue("test");
            Attr dependsCompTest = scenarioDocument.createAttribute("depends");
            dependsCompTest.setNodeValue("");
            Attr valuesCompTest = scenarioDocument.createAttribute("values");
            valuesCompTest.setNodeValue("none,microscopy,RDT");

            decisionCompTest.setAttributeNode(nameCompTest);
            decisionCompTest.setAttributeNode(dependsCompTest);
            decisionCompTest.setAttributeNode(valuesCompTest);
            decisionCompTest.setTextContent("none");
            decisionsComp.appendChild(decisionCompTest);

            complicated.appendChild(decisionsComp);

            Element treatmentsComp = (Element)scenarioDocument.createElement("treatments");
            Element treatmentComp = (Element)scenarioDocument.createElement("treatment");
            Attr nameTreatComp = scenarioDocument.createAttribute("name");
            nameTreatComp.setNodeValue("effective_treat");
            treatmentComp.setAttributeNode(nameTreatComp);

            Element scheduleComp = (Element)scenarioDocument.createElement("schedule");
            Element medicateComp = (Element)scenarioDocument.createElement("medicate");

            Attr drugComp = (Attr)scenarioDocument.createAttribute("drug");
            drugComp.setNodeValue("effective");
            Attr mgComp = (Attr)scenarioDocument.createAttribute("mg");
            mgComp.setNodeValue("1");
            Attr hourComp = (Attr)scenarioDocument.createAttribute("hour");
            hourComp.setNodeValue("0");

            medicateComp.setAttributeNode(drugComp);
            medicateComp.setAttributeNode(mgComp);
            medicateComp.setAttributeNode(hourComp);

            scheduleComp.appendChild(medicateComp);
            treatmentComp.appendChild(scheduleComp);
            treatmentsComp.appendChild(treatmentComp);

            Element treatmentCompNone = (Element)scenarioDocument.createElement("treatment");
            Attr nameTreatCompNone = (Attr)scenarioDocument.createAttribute("name");
            nameTreatCompNone.setNodeValue("none");
            treatmentCompNone.setAttributeNode(nameTreatCompNone);

            Element scheduleCompNone = (Element)scenarioDocument.createElement("schedule");
            treatmentCompNone.appendChild(scheduleCompNone);
            treatmentsComp.appendChild(treatmentCompNone);
            complicated.appendChild(treatmentsComp);

            eventScheduler.appendChild(complicated);

            Element clinicalOutcomes = (Element)scenarioDocument.createElement("ClinicalOutcomes");
            Element maxUCSeekingMemory = (Element)scenarioDocument.createElement("maxUCSeekingMemory");
            maxUCSeekingMemory.setTextContent("3");
            Element uncomplicatedCaseDuration = (Element)scenarioDocument.createElement("uncomplicatedCaseDuration");
            uncomplicatedCaseDuration.setTextContent("3");
            Element complicatedCaseDuration = (Element)scenarioDocument.createElement("complicatedCaseDuration");
            complicatedCaseDuration.setTextContent("5");
            Element complicatedRiskDuration = (Element)scenarioDocument.createElement("complicatedRiskDuration");
            complicatedRiskDuration.setTextContent("5");
            Element pImmediateUC = (Element)scenarioDocument.createElement("pImmediateUC");
            pImmediateUC.setTextContent("1");
            //Element propDeathsFirstDay = (Element)scenarioDocument.createElement("propDeathsFirstDay");
            //propDeathsFirstDay.setTextContent("0.4");

            //this communityOddsMultiplier will be removed (for schema >= 19)
            //Element communityOddsMultiplier = (Element)scenarioDocument.createElement("communityOddsMultiplier");
            //communityOddsMultiplier.setTextContent("1.5");

            clinicalOutcomes.appendChild(maxUCSeekingMemory);
            clinicalOutcomes.appendChild(uncomplicatedCaseDuration);
            clinicalOutcomes.appendChild(complicatedCaseDuration);
            clinicalOutcomes.appendChild(complicatedRiskDuration);
            clinicalOutcomes.appendChild(pImmediateUC);
            

            //clinicalOutcomes.appendChild(communityOddsMultiplier);

            eventScheduler.appendChild(clinicalOutcomes);


            Element CFR = (Element)healthSystem.getElementsByTagName("CFR").item(0);
            if(CFR!=null)
                healthSystem.insertBefore(eventScheduler, CFR);
            else
                healthSystem.appendChild(eventScheduler);
        }
        else if(valueString.equals("Ironmal"))
        {
            Element eventScheduler = (Element)scenarioDocument.createElement("EventScheduler");
            healthSystem.removeChild(immediateOutcomes);

            Element uncomplicated = (Element)scenarioDocument.createElement("uncomplicated");
            Element decisions = (Element)scenarioDocument.createElement("decisions");

            Element decisionOC = (Element)scenarioDocument.createElement("decision");

            Attr nameOC = scenarioDocument.createAttribute("name");
            nameOC.setNodeValue("official_care");
            Attr dependsOC = scenarioDocument.createAttribute("depends");
            dependsOC.setNodeValue("p");
            Attr valuesOC = scenarioDocument.createAttribute("values");
            valuesOC.setNodeValue("yes,no");

            decisionOC.setAttributeNode(nameOC);
            decisionOC.setAttributeNode(dependsOC);
            decisionOC.setAttributeNode(valuesOC);
            decisionOC.setTextContent("p(.64): yes p(.36): no");
            decisions.appendChild(decisionOC);

            Element decisionTreat = (Element)scenarioDocument.createElement("decision");

            Attr nameTreat = scenarioDocument.createAttribute("name");
            nameTreat.setNodeValue("treatment");
            Attr dependsTreat = scenarioDocument.createAttribute("depends");
            dependsTreat.setNodeValue("official_care,p");
            Attr valuesTreat = scenarioDocument.createAttribute("values");
            valuesTreat.setNodeValue("effective_treat,none");

            decisionTreat.setAttributeNode(nameTreat);
            decisionTreat.setAttributeNode(dependsTreat);
            decisionTreat.setAttributeNode(valuesTreat);
            decisionTreat.setTextContent("official_care(yes){p(.6): effective_treat p(.4): none} official_care(no): none");
            decisions.appendChild(decisionTreat);

            Element decisionTest = (Element)scenarioDocument.createElement("decision");

            Attr nameTest = scenarioDocument.createAttribute("name");
            nameTest.setNodeValue("test");
            Attr dependsTest = scenarioDocument.createAttribute("depends");
            dependsTest.setNodeValue("");
            Attr valuesTest = scenarioDocument.createAttribute("values");
            valuesTest.setNodeValue("none,microscopy,RDT");

            decisionTest.setAttributeNode(nameTest);
            decisionTest.setAttributeNode(dependsTest);
            decisionTest.setAttributeNode(valuesTest);
            decisionTest.setTextContent("none");
            decisions.appendChild(decisionTest);

            uncomplicated.appendChild(decisions);

            Element treatments = (Element)scenarioDocument.createElement("treatments");
            Element treatment = (Element)scenarioDocument.createElement("treatment");
            Attr nameTreatEl = scenarioDocument.createAttribute("name");
            nameTreatEl.setNodeValue("effective_treat");
            treatment.setAttributeNode(nameTreatEl);

            Element schedule = (Element)scenarioDocument.createElement("schedule");
            Element medicate = (Element)scenarioDocument.createElement("medicate");

            Attr drug = (Attr)scenarioDocument.createAttribute("drug");
            drug.setNodeValue("effective");
            Attr mg = (Attr)scenarioDocument.createAttribute("mg");
            mg.setNodeValue("1");
            Attr hour = (Attr)scenarioDocument.createAttribute("hour");
            hour.setNodeValue("0");

            medicate.setAttributeNode(drug);
            medicate.setAttributeNode(mg);
            medicate.setAttributeNode(hour);

            schedule.appendChild(medicate);
            treatment.appendChild(schedule);
            treatments.appendChild(treatment);

            Element treatmentNone = (Element)scenarioDocument.createElement("treatment");
            Attr nameTreatNone = (Attr)scenarioDocument.createAttribute("name");
            nameTreatNone.setNodeValue("none");
            treatmentNone.setAttributeNode(nameTreatNone);

            Element scheduleNone = (Element)scenarioDocument.createElement("schedule");
            treatmentNone.appendChild(scheduleNone);
            treatments.appendChild(treatmentNone);
            uncomplicated.appendChild(treatments);

            eventScheduler.appendChild(uncomplicated);

            Element complicated = (Element)scenarioDocument.createElement("complicated");
            Element decisionsComp = (Element)scenarioDocument.createElement("decisions");

            Element decisionOCComp = (Element)scenarioDocument.createElement("decision");

            Attr nameOCComp = scenarioDocument.createAttribute("name");
            nameOCComp.setNodeValue("official_care");
            Attr dependsOCComp = scenarioDocument.createAttribute("depends");
            dependsOCComp.setNodeValue("p");
            Attr valuesOCComp = scenarioDocument.createAttribute("values");
            valuesOCComp.setNodeValue("yes,no");

            decisionOCComp.setAttributeNode(nameOCComp);
            decisionOCComp.setAttributeNode(dependsOCComp);
            decisionOCComp.setAttributeNode(valuesOCComp);
            decisionOCComp.setTextContent("p(.48): yes p(.52): no");
            decisionsComp.appendChild(decisionOCComp);

            Element decisionCompTreat = (Element)scenarioDocument.createElement("decision");

            Attr nameCompTreat = scenarioDocument.createAttribute("name");
            nameCompTreat.setNodeValue("treatment");
            Attr dependsCompTreat = scenarioDocument.createAttribute("depends");
            dependsCompTreat.setNodeValue("official_care,p");
            Attr valuesCompTreat = scenarioDocument.createAttribute("values");
            valuesCompTreat.setNodeValue("effective_treat,none");

            decisionCompTreat.setAttributeNode(nameCompTreat);
            decisionCompTreat.setAttributeNode(dependsCompTreat);
            decisionCompTreat.setAttributeNode(valuesCompTreat);
            decisionCompTreat.setTextContent("official_care(yes){p(.6): effective_treat p(.4): none} official_care(no): none");
            decisionsComp.appendChild(decisionCompTreat);

            Element decisionCompHosp = (Element)scenarioDocument.createElement("decision");

            Attr nameCompHosp = scenarioDocument.createAttribute("name");
            nameCompHosp.setNodeValue("hospitalisation");
            Attr dependsCompHosp = scenarioDocument.createAttribute("depends");
            dependsCompHosp.setNodeValue("official_care");
            Attr valuesCompHosp = scenarioDocument.createAttribute("values");
            valuesCompHosp.setNodeValue("none,delayed,immediate");

            decisionCompHosp.setAttributeNode(nameCompHosp);
            decisionCompHosp.setAttributeNode(dependsCompHosp);
            decisionCompHosp.setAttributeNode(valuesCompHosp);
            decisionCompHosp.setTextContent("official_care(yes): immediate official_care(no): none");
            decisionsComp.appendChild(decisionCompHosp);

            Element decisionCompTest = (Element)scenarioDocument.createElement("decision");

            Attr nameCompTest = scenarioDocument.createAttribute("name");
            nameCompTest.setNodeValue("test");
            Attr dependsCompTest = scenarioDocument.createAttribute("depends");
            dependsCompTest.setNodeValue("");
            Attr valuesCompTest = scenarioDocument.createAttribute("values");
            valuesCompTest.setNodeValue("none,microscopy,RDT");

            decisionCompTest.setAttributeNode(nameCompTest);
            decisionCompTest.setAttributeNode(dependsCompTest);
            decisionCompTest.setAttributeNode(valuesCompTest);
            decisionCompTest.setTextContent("none");
            decisionsComp.appendChild(decisionCompTest);

            complicated.appendChild(decisionsComp);

            Element treatmentsComp = (Element)scenarioDocument.createElement("treatments");
            Element treatmentComp = (Element)scenarioDocument.createElement("treatment");
            Attr nameTreatComp = scenarioDocument.createAttribute("name");
            nameTreatComp.setNodeValue("effective_treat");
            treatmentComp.setAttributeNode(nameTreatComp);

            Element scheduleComp = (Element)scenarioDocument.createElement("schedule");
            Element medicateComp = (Element)scenarioDocument.createElement("medicate");

            Attr drugComp = (Attr)scenarioDocument.createAttribute("drug");
            drugComp.setNodeValue("effective");
            Attr mgComp = (Attr)scenarioDocument.createAttribute("mg");
            mgComp.setNodeValue("1");
            Attr hourComp = (Attr)scenarioDocument.createAttribute("hour");
            hourComp.setNodeValue("0");

            medicateComp.setAttributeNode(drugComp);
            medicateComp.setAttributeNode(mgComp);
            medicateComp.setAttributeNode(hourComp);

            scheduleComp.appendChild(medicateComp);
            treatmentComp.appendChild(scheduleComp);
            treatmentsComp.appendChild(treatmentComp);

            Element treatmentCompNone = (Element)scenarioDocument.createElement("treatment");
            Attr nameTreatCompNone = (Attr)scenarioDocument.createAttribute("name");
            nameTreatCompNone.setNodeValue("none");
            treatmentCompNone.setAttributeNode(nameTreatCompNone);

            Element scheduleCompNone = (Element)scenarioDocument.createElement("schedule");
            treatmentCompNone.appendChild(scheduleCompNone);
            treatmentsComp.appendChild(treatmentCompNone);
            complicated.appendChild(treatmentsComp);

            eventScheduler.appendChild(complicated);

            Element clinicalOutcomes = (Element)scenarioDocument.createElement("ClinicalOutcomes");
            Element maxUCSeekingMemory = (Element)scenarioDocument.createElement("maxUCSeekingMemory");
            maxUCSeekingMemory.setTextContent("3");
            Element uncomplicatedCaseDuration = (Element)scenarioDocument.createElement("uncomplicatedCaseDuration");
            uncomplicatedCaseDuration.setTextContent("3");
            Element complicatedCaseDuration = (Element)scenarioDocument.createElement("complicatedCaseDuration");
            complicatedCaseDuration.setTextContent("5");
            Element complicatedRiskDuration = (Element)scenarioDocument.createElement("complicatedRiskDuration");
            complicatedRiskDuration.setTextContent("5");
            Element pImmediateUC = (Element)scenarioDocument.createElement("pImmediateUC");
            pImmediateUC.setTextContent("1");
            //Element propDeathsFirstDay = (Element)scenarioDocument.createElement("propDeathsFirstDay");
            //propDeathsFirstDay.setTextContent("0.4");
            
            //this communityOddsMultiplier will be removed (for schema >= 19)
            //Element communityOddsMultiplier = (Element)scenarioDocument.createElement("communityOddsMultiplier");
            //communityOddsMultiplier.setTextContent("1.5");

            clinicalOutcomes.appendChild(maxUCSeekingMemory);
            clinicalOutcomes.appendChild(uncomplicatedCaseDuration);
            clinicalOutcomes.appendChild(complicatedCaseDuration);
            clinicalOutcomes.appendChild(complicatedRiskDuration);
            clinicalOutcomes.appendChild(pImmediateUC);
            

            //clinicalOutcomes.appendChild(communityOddsMultiplier);

            eventScheduler.appendChild(clinicalOutcomes);


            Element CFR = (Element)healthSystem.getElementsByTagName("CFR").item(0);
            if(CFR!=null)
                healthSystem.insertBefore(eventScheduler, CFR);
            else
                healthSystem.appendChild(eventScheduler);
        }
        else return false;

        // creating drugdescription element
        Element drugDescription = scenarioDocument.createElement("drugDescription");

        Element drug = scenarioDocument.createElement("drug");
        Attr abbrev = scenarioDocument.createAttribute("abbrev");
        abbrev.setNodeValue("effective");
        drug.setAttributeNode(abbrev);

        Element pd = scenarioDocument.createElement("PD");

        Element allele = scenarioDocument.createElement("allele");
        Attr nameAllele = scenarioDocument.createAttribute("name");
        nameAllele.setNodeValue("sensitive");
        allele.setAttributeNode(nameAllele);

        Element initial_frequency = scenarioDocument.createElement("initial_frequency");
        initial_frequency.setTextContent("1");

        Element max_killing_rate = scenarioDocument.createElement("max_killing_rate");
        max_killing_rate.setTextContent("1e7");

        Element ic50 = scenarioDocument.createElement("IC50");
        ic50.setTextContent("1");

        Element slope = scenarioDocument.createElement("slope");
        slope.setTextContent("1");

        allele.appendChild(initial_frequency);
        allele.appendChild(max_killing_rate);
        allele.appendChild(ic50);
        allele.appendChild(slope);

        pd.appendChild(allele);
        drug.appendChild(pd);


        Element pk = scenarioDocument.createElement("PK");

        Element negligible_concentration = scenarioDocument.createElement("negligible_concentration");
        negligible_concentration.setTextContent("1e-5");

        Element half_life = scenarioDocument.createElement("half_life");
        half_life.setTextContent("0.00069");

        Element vol_dist = scenarioDocument.createElement("vol_dist");
        vol_dist.setTextContent("0.01667");

        pk.appendChild(negligible_concentration);
        pk.appendChild(half_life);
        pk.appendChild(vol_dist);
        drug.appendChild(pk);

        drugDescription.appendChild(drug);

        scenarioElement.insertBefore(drugDescription, modelElement);

    	return true;
    }
    
    private void setMosqsNewAttributes(int mosqType, Element mosq,
            NodeList nonHumanHosts) {

        Attr humanBloodIndex = scenarioDocument
                .createAttribute("mosqHumanBloodIndex");
        Attr proportionLaidEggsSameDay = scenarioDocument
                .createAttribute("mosqLaidEggsSameDayProportion");
        Attr PSurvivalFeedingCycle = scenarioDocument
                .createAttribute("mosqSurvivalFeedingCycleProbability");

        proportionLaidEggsSameDay
                .setNodeValue(Double
                        .toString(SchemaTranslator.ProporitionsLaidEggsSameDay[mosqType]));
        PSurvivalFeedingCycle.setNodeValue(Double
                .toString(SchemaTranslator.PsSurvivalFeedingCycle[mosqType]));

        if (nonHumanHosts == null || nonHumanHosts.getLength() == 0)
            humanBloodIndex.setNodeValue(Double
                    .toString(SchemaTranslator.HumanBloodIndex_NONNHS));

        else if (nonHumanHosts.getLength() == 1) {
            humanBloodIndex.setNodeValue(Double
                    .toString(SchemaTranslator.HumanBloodIndexes[mosqType]));

            Element nhh = (Element) nonHumanHosts.item(0);

            Attr relativeEntoAvailability = scenarioDocument
                    .createAttribute("mosqRelativeEntoAvailability");
            relativeEntoAvailability.setNodeValue(Double
                    .toString(SchemaTranslator.Standard_RELATIVE_ENTO_AV));

            nhh.setAttributeNode(relativeEntoAvailability);
            nhh.removeAttribute("mosqEntoAvailability");
        } else {

            humanBloodIndex.setNodeValue(Double
                    .toString(SchemaTranslator.HumanBloodIndexes[mosqType]));
            System.err
                    .println("There are more than 1 non human hosts types in these scenario. Please edit the relative Ento availabilities for each type of non human host by hand.");
        }

        mosq.setAttributeNode(humanBloodIndex);
        mosq.setAttributeNode(proportionLaidEggsSameDay);
        mosq.setAttributeNode(PSurvivalFeedingCycle);

        mosq.removeAttribute("mosqEntoAvailability");
        mosq.removeAttribute("mosqSeekingDeathRate");
    }

    private void visitAllFiles(File file, File outDir) throws Exception {
        if (file.isDirectory()) {
            String[] children = file.list();
            for (int i = 0; i < children.length; i++) {
                try {
                    visitAllFiles(new File(file, children[i]), outDir);
                } catch (Exception exc) {
                    System.out
                            .println("Error translating " + file + ": " + exc);
                    exc.printStackTrace();
                }
            }
        } else {
            if (file.getName().endsWith(".xml")) {
                System.out.println(file.getAbsolutePath());
                translateFile(file, outDir);
            }
        }
    }

    public static void main(String[] args) {
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--required_version")) {
                _required_version = Integer.parseInt(args[++i]);
            } else if (args[i].equals("--oneDayTimesteps")) {
            	doODTTranslation = true;
                doValidation = false;
                System.out.println("You have chosen the --oneDayTimesteps option, this option is only intended for the fitting scenarii or scenarii using no intervention and/or no ");
            } else if (args[i].equals("--no-validation")) {
                doValidation = false;
            } else if (args[i].equals("--no-translation")) {
                doTranslation = false;
            } else if (args[i].equals("--update-db")) {
                doDBUpdate = true;
            } else if (args[i].equals("--maxDensCorrection")) {
                String arg = args[++i];
                if (arg.equalsIgnoreCase("true")) {
                    maxDensBug = BugCorrectionBehaviour.correct;
                } else if (arg.equalsIgnoreCase("false")) {
                    maxDensBug = BugCorrectionBehaviour.dontCorrect;
                } else {
                    System.err
                            .println("--maxDensCorrection: expected true or false");
                    System.exit(2);
                }
            } else {
                printUsage();
            }
            System.out.println(args[i]);
        }
        if (_required_version == 1) {
            System.out
                    .println("Target version 1 is not supported (see comment for translate0To1).");
            System.exit(1);
        }

        SchemaTranslator st = new SchemaTranslator();
        try {
            if (doDBUpdate) {
                st.updateDB();
            } else {
                File scenarios = new File("scenarios");
                if (!scenarios.isDirectory())
                    scenarios.mkdir();
                System.out
                        .println("Put XMLs to be translated into the \"scenarios\" directory");
                File outDir = new File("translatedScenarios");
                if (!outDir.isDirectory())
                    outDir.mkdir();
                st.visitAllFiles(scenarios, outDir);
            }

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    private static void printUsage() {
        System.out
                .println("Usage: schemaTranslator [options]:\n"
                        + "--required_version VERSION\tThe version number to update the document(s) to. Default: CURRENT_VERSION"
                        + "--no-validation\t\tDon't validate the result"
                        + "--no-translation\t\tDon't write out the translated result (but still translate internally for validation)"
                        + "--update-db\t\tUpdate DB entries instead of files"
                        + "--maxDensCorrection BOOL\tUpdate 12->13 requires this sometimes: set true to include bug fix, false to explicitly exclude it.");
        System.exit(1);
    }
}
