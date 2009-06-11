import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.lang.reflect.Method;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Result;
import javax.xml.transform.Source;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
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
import org.xml.sax.SAXException;

public class SchemaTranslator {

	DocumentBuilder _builder;
	Document scenarioDocument;
	Element scenarioElement;

	static final int CURRENT_VERSION = 5;

	private static int _required_version = CURRENT_VERSION;
	private static boolean doValidation = true;
	private static boolean doTranslation = true;

	public SchemaTranslator() {
		try {
			_builder = DocumentBuilderFactory.newInstance()
					.newDocumentBuilder();
		} catch (ParserConfigurationException e) {
			e.printStackTrace();
		}
	}

	private void validate(Document scenarioDocument, String schemaFileName) throws SAXException,
			IOException {
		Document forValidation = (Document) scenarioDocument.cloneNode(true);
		Element scenarioElement = forValidation.getDocumentElement();
		scenarioElement.removeAttribute("xsi:noNamespaceSchemaLocation");

		// This is set by generateRun
		scenarioElement.setAttribute("assimMode", "0");
		// This is set by the work generator
		scenarioElement.setAttribute("wuID", "123");

		if (scenarioElement.getLastChild().getNodeValue() == "@parameters@") {
			scenarioElement.getLastChild().setNodeValue("");
			// Add a dummy parameter
			Element parameters = forValidation.createElement("parameters");

			parameters.setAttribute("latentp", "0");
			parameters.setAttribute("delta", "0");
			parameters.setAttribute("eipDuration", "0");
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
		File xsdFile = new File("../../test/original/"+schemaFileName);
		if (xsdFile == null || !xsdFile.isFile()) {
			System.out.println("Unable to find scenario.xsd file; not validating.");
			return;
		}
		SchemaFactory factory = SchemaFactory
				.newInstance("http://www.w3.org/2001/XMLSchema");
		//NOTE: this may throw a java.lang.NullPointerException when run within eclipse
		// run from a command-line instead (java SchemaTranslator).
		Schema schema = factory.newSchema(xsdFile);
		Validator validator = schema.newValidator();

		Source source = new DOMSource(forValidation);
		validator.validate(source);
	}

	private void translate(File documentFile, File outDir) throws Exception {
		scenarioDocument = _builder.parse(documentFile);
		scenarioElement = scenarioDocument.getDocumentElement();
		// 0 if no current version (getAttribute returns ""):
		int schemaVersion = Integer.parseInt ("0"+scenarioElement.getAttribute("schemaVersion"));
		String schemaFileName = "scenario"+schemaVersion+".xsd";
		Class<SchemaTranslator> cls = (Class<SchemaTranslator>)this.getClass();
		while (schemaVersion < _required_version) {
			++schemaVersion;
			schemaFileName = "scenario_"+schemaVersion+".xsd";
			scenarioElement.setAttribute("schemaVersion", Integer.toString(schemaVersion));
			scenarioElement.setAttribute("xsi:noNamespaceSchemaLocation",
			schemaFileName);
			
			String translateMeth = "translate"+(schemaVersion-1)+"To"+schemaVersion;
			Method method = cls.getMethod(translateMeth, new Class[] {});
			method.invoke(this, new Object[] {});
		}
		if (doTranslation) {
			File outFile = new File(outDir, documentFile.getName());
			outFile.createNewFile();
			OutputStream os = new FileOutputStream(outFile);
			/*Result result = new StreamResult(os);
			// Write the DOM document to the file
			Transformer xformer = TransformerFactory.newInstance().newTransformer();
			xformer.setOutputProperty(OutputKeys.ENCODING, "UTF-8");
			xformer.setOutputProperty(OutputKeys.METHOD, "xml");
			// This adds more indentation/new-lines where there's already spacing :-( 
			//xformer.setOutputProperty(OutputKeys.INDENT, "yes");
			xformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "2");
			xformer.transform(new DOMSource(scenarioDocument), result);*/
	        OutputFormat format = new OutputFormat(doc);
	         format.setLineWidth(65);
	         format.setIndenting(true);
	         format.setIndent(2);
	         XMLSerializer serializer = new XMLSerializer(out, format);
	         serializer.serialize(doc);
	 		}
		if (doValidation) validate(scenarioDocument, schemaFileName);
	}

	/// Exactly what version 1 is has been forgotten; it's merged into 2.
	public void translate0To1() {
	}
	public void translate1To2() {
		scenarioElement.setAttribute("xmlns:xsi",
				"http://www.w3.org/2001/XMLSchema-instance");
		scenarioElement.setAttribute("xsi:noNamespaceSchemaLocation",
				"scenario.xsd");
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
	}

	public void translate2To3() {
	}
	
	/* EntoData now has either a nonVector or a vector element; EIRDaily and
	 * anopheles lists have moved to one of these.
	 * Some unwanted entomological parameters have been removed, many have been
	 * added, and eipDuration has been moved from param. */
	public void translate3To4() {
		Element params = (Element)scenarioElement.getElementsByTagName("parameters").item(0);
		Attr eip = params.getAttributeNode("eipDuration");
		params.removeAttributeNode(eip);
		// Attribute added to nonVector, if used, below.
		
		Element elt = (Element) scenarioElement
		.getElementsByTagName("entoData").item(0);
		if (elt != null) {
			NodeList list = elt.getElementsByTagName("EIRDaily");
			if (list.getLength() != 0) {
				Element nonVector = scenarioDocument.createElement("nonVector");
				elt.appendChild(nonVector);
				// NOTE: results don't seem to be right if iterating forwards instead of backwards
				for (int i = list.getLength()-1; i >= 0; --i) {
					//System.out.println (list.getLength());
					Node eir = list.item(i);
					// add the element at the new location, remove from old:
					// unfortunately messes up the whitespace
					nonVector.insertBefore(eir, nonVector.getFirstChild());
				}
				nonVector.setAttributeNode(eip);
			}
			list = elt.getElementsByTagName("anopheles");
			if (list.getLength() != 0) {
				Element vector = scenarioDocument.createElement("vector");
				elt.appendChild(vector);
				for (int i = list.getLength()-1; i >= 0; --i) {
					Node anoph = list.item(i);
					((Element)anoph).removeAttribute("useNv0Guess");
					vector.insertBefore(anoph, vector.getFirstChild());
				}
			}
			
			elt.removeAttribute("inputType");
		}
		
		elt = (Element) scenarioElement.getElementsByTagName("interventions").item(0);
		if (elt != null) {
			Element timed = (Element) elt.getElementsByTagName("timed").item(0);
			if (timed != null) {
				Element interv = (Element) timed.getElementsByTagName("intervention").item(0);
				if (interv != null) {
					Element cEIR = (Element) interv.getElementsByTagName("changeEIR").item(0);
					if (cEIR != null) {
						cEIR.removeAttribute("inputType");
						cEIR.removeAttribute("name");	// part of EntoData not NonVector
						cEIR.setAttributeNode((Attr)eip.cloneNode(true));
					}
				}
			}
		}
	}
	
	// modelVersion flags 1<<2, 1<<4 or 1<<5 have changed
	public void translate4To5() throws Exception {
		int ver = Integer.parseInt(scenarioElement.getAttribute("modelVersion"));
		if ((ver & 0x68) != 0)	// modelVersion with flags 1<<2, 1<<4 or 1<<5
			throw new Exception ("Scenario uses InfectionIncidence model flags; these have changed. Please update!");
		/* Note : We don't translate directly because
		 *  a) Very few scenarios already do this
		 *  b) It's probably better that people know this change is required.
		 * However, there is a direct translation in 2/3 of cases:
		 *  1<<2 translates to (1<<2 | 1<<5),
		 *  1<<4 translates to (1<<4 | 1<<5),
		 *  1<<5 had a bug fixed, but roughly translates to 1<<4.
		 */
	}

	private void visitAllFiles(File file, File outDir) throws Exception {
		if (file.isDirectory()) {
			String[] children = file.list();
			for (int i = 0; i < children.length; i++) {
				try {
					visitAllFiles(new File(file, children[i]), outDir);
				} catch (Exception exc) {
					System.out.println ("Error translating "+file+": "+exc);
				}
			}
		} else {
			if (file.getName().endsWith(".xml")) {
				System.out.println(file.getAbsolutePath());
				translate(file, outDir);
			}
		}
	}

	public static void main(String[] args) {
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("--required_version")) {
				_required_version = Integer.parseInt(args[++i]);
			} else if (args[i].equals("--no-validation")) {
				doValidation = false;
			} else if (args[i].equals("--no-translation")) {
				doTranslation = false;
			} else {
				printUsage();
			}
			System.out.println(args[i]);
		}
		if (_required_version == 1) {
			System.out.println("Target version 1 is not supported (see comment for translate0To1).");
			System.exit(1);
		}

		SchemaTranslator st = new SchemaTranslator();
		try {
			File scenarios = new File("scenarios");
			if (!scenarios.isDirectory())
				scenarios.mkdir();
			File outDir = new File("translatedScenarios");
			if (!outDir.isDirectory())
				outDir.mkdir();
			st.visitAllFiles(scenarios, outDir);
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	private static void printUsage() {
		System.out.println("Usage: schemaTranslator [options]:\n"
				+ "--required_version VERSION\tThe version number to update the document(s) to. Default: CURRENT_VERSION"
				+ "--no-validation\t\tDon't validate the result"
				+ "--no-translation\t\tDon't write out the translated result (but still translate internally for validation)"
				);
		System.exit(1);
	}
}
