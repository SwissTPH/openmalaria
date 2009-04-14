import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.Result;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

public class SchemaTranslator {

    DocumentBuilder _builder;
    static final int CURRENT_VERSION = 3;
    private int _required_version;

    public SchemaTranslator(int required_version) {
        _required_version = required_version;
        try {
            _builder = DocumentBuilderFactory.newInstance()
                    .newDocumentBuilder();
        } catch (ParserConfigurationException e) {
            e.printStackTrace();
        }
    }

    private void translate(File documentFile) throws Exception {
        Document scenarioDocument = _builder.parse(documentFile);
        Element scenarioElement = scenarioDocument.getDocumentElement();
        String schemaVersion = scenarioElement.getAttribute("schemaVersion");
        if (schemaVersion == "") {
            schemaVersion = "0";
        }
        switch (Integer.parseInt(schemaVersion)) {
        case 0:
        case 1:
            translate0To2(scenarioElement);
            if (_required_version == Integer.parseInt(scenarioElement
                    .getAttribute("schemaVersion"))) {
                break;
            }
        case 2:
            translate2To3(scenarioElement);
            if (_required_version == Integer.parseInt(scenarioElement
                    .getAttribute("schemaVersion"))) {
                break;
            }
        case 3:
            break;
        default:
            System.out.println("Not a recognized schemaVersion");
        }
        File outFile = new File("translatedScenarios/" + documentFile.getName());
        OutputStream os = new FileOutputStream(outFile);
        Result result = new StreamResult(os);
        // Write the DOM document to the file
        Transformer xformer = TransformerFactory.newInstance().newTransformer();
        xformer.transform(new DOMSource(scenarioDocument), result);

    }

    private void translate0To2(Element scenarioElement) {
        scenarioElement.setAttribute("xmlns:xsi",
                "http://www.w3.org/2001/XMLSchema-instance");
        scenarioElement.setAttribute("xsi:noNamespaceSchemaLocation",
                "scenario.xsd");
        scenarioElement.setAttribute("schemaVersion", "2");
        if (((Element) scenarioElement.getElementsByTagName("entoData").item(0))
                .hasAttribute("firstDay")) {
            System.out.println("Warning: Removed firstDay attribute");
            ((Element) scenarioElement.getElementsByTagName("entoData").item(0))
                    .removeAttribute("firstDay");
        }

        if (((Element) scenarioElement.getElementsByTagName("changeEIR")
                .item(0)) != null) {
            if (((Element) scenarioElement.getElementsByTagName("changeEIR")
                    .item(0)).hasAttribute("firstDay")) {
                System.out.println("Warning: Removed firstDay attribute");
                ((Element) scenarioElement.getElementsByTagName("changeEIR")
                        .item(0)).removeAttribute("firstDay");
            }
        }
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
    }

    private void translate2To3(Element scenarioElement) {
        scenarioElement.setAttribute("schemaVersion", "3");
    }

    private void visitAllFiles(File file) throws Exception {
        if (file.isDirectory()) {
            String[] children = file.list();
            for (int i = 0; i < children.length; i++) {
                visitAllFiles(new File(file, children[i]));
            }
        } else {
            if (file.getName().endsWith(".xml")) {
                System.out.println(file.getAbsolutePath());
                translate(file);
            }
        }
    }

    public static void main(String[] args) {

        int required_version = CURRENT_VERSION;

        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--required_version")) {
                required_version = Integer.parseInt(args[++i]);
            } else {
                printUsage();
            }
            System.out.println(args[i]);
        }

        SchemaTranslator st = new SchemaTranslator(required_version);
        try {
            st.visitAllFiles(new File("scenarios"));
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    private static void printUsage() {
        System.out.println("Usage: Supply the following cmdl options:\n"
                + "--required_version (optional): The version number to update"
                + "the document(s) to. Default: CURRENT_VERSION");
        System.exit(1);
    }
}
