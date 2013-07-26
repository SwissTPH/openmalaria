package SchemaTranslator

import java.io.File
import javax.xml.parsers.DocumentBuilderFactory
import org.w3c.dom.Document
import org.xml.sax.ErrorHandler
import org.xml.sax.InputSource
import org.w3c.dom.Element
import org.w3c.dom.Node
import org.w3c.dom.NodeList
import kotlin.dom.*
import java.util.ArrayList
import javax.xml.transform.TransformerFactory
import javax.xml.transform.Result
import javax.xml.transform.OutputKeys
import javax.xml.transform.dom.DOMSource
import javax.xml.validation.SchemaFactory
import javax.xml.validation.Schema
import javax.xml.validation.Validator
import javax.xml.transform.Source
import org.xml.sax.SAXParseException
import java.lang.reflect.Method
import java.lang.reflect.InvocationTargetException
import java.util.TreeSet

// ———  part 1: define utility classes (enums, simple containers)  ———

class SetupException(msg: String) : Exception(msg) {}

class DocumentException(msg: String) : Exception(msg) {}

open class STErrorHandler() : ErrorHandler {
    public override fun fatalError(p0: SAXParseException?) : Unit {
        System.err.println("Error: " + p0!!.toString())
        throw p0
    }
    public override fun error(p0: SAXParseException?) : Unit {
        System.err.println("Error: " + p0!!.toString())
    }
    public override fun warning(p0: SAXParseException?) : Unit {
        System.err.println("Error: " + p0!!.toString())
    }
}

enum class SchemaName {
    VERSIONED; NO_SUFFIX; CURRENT
}
enum class BugCorrectionBehaviour {
    NONE; CORRECT; DONT_CORRECT
}
enum class IptiSpBehaviour {
    NONE; ASSUME_INTENDED; ASSUME_UNINTENDED
}
enum class IptiReportOnlyAtRiskBehaviour {
    NONE; ON; OFF
}
enum class ITN29ParameterTranslation {
    NONE; REPLACE; MANUAL  // note: could add option to approximate old behaviour
}

class Options {
    var outputFolder = File("translatedScenarios")
    var inputFolder = File("scenarios")
    var schemaFolder = File("../../schema/")

    val latestVersion = 32
    var targetVersion = latestVersion

    var doValidation = true
    var doTranslation = true
    var doODTTranslation = false
    var doDBUpdate = false

    var latestSchema = SchemaName.VERSIONED

    var maxDensBug = BugCorrectionBehaviour.NONE

    var iptiSpOption = IptiSpBehaviour.NONE

    var iptiROAR = IptiReportOnlyAtRiskBehaviour.NONE

    var ITN29Translation = ITN29ParameterTranslation.NONE;
}

val transformerFactory = TransformerFactory.newInstance()!!

// ———  part 2: translation class and utility functions  ———

fun validate(scenarioDocument: Document, schemaFile: File, options: Options) : Unit {
    val forValidation: Document = scenarioDocument.cloneNode(true) as Document
    val scenarioElement: Element = forValidation.getDocumentElement()!!
    scenarioElement.removeAttribute("xsi:noNamespaceSchemaLocation")
    if (options.targetVersion <= 23){
        scenarioElement.setAttribute("assimMode", "0")
    }

    scenarioElement.setAttribute("wuID", "123")
    val t_parameters: Element? = scenarioElement.getElementsByTagName("parameters").item(0) as Element?
    if (t_parameters != null && t_parameters.getNodeValue() != null &&
    (t_parameters.getNodeValue()!!.contains("@parameters@")))
    {
        t_parameters.getLastChild()!!.setNodeValue("")
        val parameters: Element = forValidation.createElement("parameters")!!
        parameters.setAttribute("latentp", "0")
        parameters.setAttribute("delta", "0")
        parameters.setAttribute("interval", "0")
        parameters.setAttribute("iseed", "0")
        val parameter: Element = forValidation.createElement("parameter")!!
        parameter.setAttribute("value", "0")
        parameter.setAttribute("name", "0")
        parameter.setAttribute("number", "0")
        parameters.appendChild(parameter)
        scenarioElement.appendChild(parameters)
    }

    val factory: SchemaFactory = SchemaFactory.newInstance("http://www.w3.org/2001/XMLSchema")
    val schema: Schema = factory.newSchema(schemaFile)!!
    val validator: Validator = schema.newValidator()!!
    validator.setErrorHandler(STErrorHandler())
    val source: Source = DOMSource(forValidation)
    validator.validate(source)
}

abstract class Translator(input: InputSource, options: Options) {
    protected val options: Options = options
    val builder = DocumentBuilderFactory.newInstance()?.newDocumentBuilder()!!
    protected var scenarioDocument: Document = builder.parse(input)!!
    protected var scenarioElement: Element = scenarioDocument.getDocumentElement()!!

    fun writeTo(result: Result){
        // Write the DOM document to the file
        val xformer = transformerFactory.newTransformer()!!
        xformer.setOutputProperty(OutputKeys.ENCODING, "UTF-8")
        xformer.setOutputProperty(OutputKeys.METHOD, "xml")

        // Reformat. First, remove all whitespace nodes, then insert new whitespace.
        stripWhitespace(scenarioElement, org.w3c.dom.Node.TEXT_NODE, "#text")
        scenarioDocument.normalize()
        // Then add new indentation/new-lines:
        xformer.setOutputProperty(OutputKeys.INDENT, "yes")
        xformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "2")

        xformer.transform(DOMSource(scenarioDocument), result)
    }
    // Helper function to strip old white-space.
    // Source: http://forums.sun.com/thread.jspa?threadID=5201482
    fun stripWhitespace(start: Node, nodeType: Short, name: String) {
        fun visit(node: Node){
            for (child in node.getChildNodes().toList()){
                if (child.getNodeType() == nodeType &&
                child.getNodeName()!!.trim() == name &&
                child.getNodeValue()!!.trim() == "")
                {
                    child.getParentNode()!!.removeChild(child)
                    // child was removed so list invalid; easiest is to start again:
                    visit(node)
                    break;
                } else {
                    visit(child)
                }
            }
        }
        visit(start)
    }

    fun getChildNodes(node: Node, name: String) : ArrayList<Node> {
        var children: NodeList = node.getChildNodes()
        var l: Int = children.getLength()
        val r = ArrayList<Node>()
        for (i in 0..l - 1) {
            val item = children.item(i)!!
            if (name == item.getNodeName())
                r.add(item)
        }
        return r
    }
    fun getChildElements(node: Node, name: String) : ArrayList<Element> {
        var children: NodeList = node.getChildNodes()
        var l: Int = children.getLength()
        val r = ArrayList<Element>()
        for (i in 0..l - 1) {
            val item = children.item(i)!!
            if (item is Element && name == item.getNodeName())
                r.add(item)
        }
        return r
    }
    fun getChildElement(node: Node, name: String): Element {
        var elts : List<Node> = getChildNodes(node, name)
        if (elts.size > 1)
            throw DocumentException("Expected ${node.getNodeName()} not to have more than one sub-element with name ${name}")
        if (elts.size == 1) return elts.get(0) as Element
        throw DocumentException( "Node ${node.getNodeName()} does not have required child ${name}" )
    }
    fun getChildElementOpt(node: Node, name: String): Element? {
        var elts : List<Node> = getChildNodes(node, name)
        if (elts.size > 1)
            throw DocumentException("Expected ${node.getNodeName()} not to have more than one sub-element with name ${name}")
        return if (elts.size == 1) elts.get(0) as Element else null
    }
    fun getOrCreateSubElt(parent: Element, name: String): Element {
        val child: Element? = getChildElement(parent, name)
        if (child == null){
            val child2: Element = scenarioDocument.createElement(name)!!
            parent.appendChild(child2)
            return child2
        }else
            return child
    }
    fun usesOption(name: String): Boolean {
        var opts : Element = scenarioElement.getElementsByTagName("ModelOptions").item(0) as Element
        for (n: Node in getChildNodes(opts, "option")){
            var e: Element = n as Element
            if (e.getAttribute("name") == name){
                return java.lang.Boolean.parseBoolean(e.getAttribute("value"))
            }
        }
        return name == "MAX_DENS_CORRECTION"
    }

    private fun genSchemaName(schemaVersion : Int): String {
        return when(options.latestSchema){
            SchemaName.NO_SUFFIX -> "scenario.xsd"
            SchemaName.CURRENT -> "scenario_current.xsd"
            else -> "scenario_" + schemaVersion + ".xsd"
        }
    }

    /** Translate the current document to the target version, potentially
     * also translating to use the one-day timestep.
     *
     * Throws a DocumentException on failure due to errors in the document (or
     * potentially the translation code).
     */
    fun translateAndValidate() : String {
        scenarioElement = scenarioDocument.getDocumentElement()!!
        var schemaVersion : Int = Integer.parseInt("0" + scenarioElement.getAttribute("schemaVersion"))
        var schemaFileName : String = genSchemaName(schemaVersion)
        val cls : Class<in Translator> = this.javaClass

        while (schemaVersion < options.targetVersion){
            ++schemaVersion
            schemaFileName = genSchemaName(schemaVersion)
            scenarioElement.setAttribute("schemaVersion", Integer.toString(schemaVersion))
            scenarioElement.setAttribute("xsi:noNamespaceSchemaLocation", schemaFileName)

            val translateMeth : String = "translate" + (schemaVersion - 1) + "To" + schemaVersion
            val method : Method? = cls.getMethod(translateMeth)
            if (method == null)
                throw Exception("Method " + translateMeth + " not found")
            try{
                val result = method.invoke(this) // invoke, may throw
                when (result){
                    // for backwards compatibility:
                    is Boolean -> if (!result) throw DocumentException("Translation failed (no message)")
                    // preferred:
                    null /* returned by invoke when return type is void */ -> {/*do nothing*/}
                    else -> throw Exception("Unexpected result value while updating to version ${schemaVersion}:  \"${result}\"")
                }
            }catch(e: InvocationTargetException){
                // The called method threw an exception and invoke() wrapped
                // it. Unwrap it and rethrow:
                throw e.getCause()!!
            }
        }
        if (schemaVersion == 18 && options.doODTTranslation){
            oDTTranslation()
        }

        if (options.doValidation)
            validate(scenarioDocument, File(options.schemaFolder, schemaFileName), options);

        return schemaFileName
    }

    /** One-day timestep translation.
     * 
     * Use throwDocumentException to throw since we can't add a throws clause
     * in Kotlin.
     */
    protected abstract fun oDTTranslation(): Unit

    /** Throw an exception, but pretend (for the purposes of the Java
     * compiler's exception checking) that we don't.
     * 
     * If calling a no-return function confuses the compiler too much, follow
     * the call with "return;". */
    protected fun throwUnchecked(e: Throwable): Unit{
        throw e
    }
    /** Convenient version of throwUnchecked. */
    protected fun throwDocumentException(msg: String): Unit{
        throw DocumentException(msg)
    }
}


// ———  part 3: main update routines  ———

/** Extension to hold translation functions written in Kotlin. */
abstract class TranslatorKotlin(input: InputSource, options: Options) : Translator(input,options){
    public fun translate31To32(){
        val interventions = getChildElement(scenarioElement, "interventions")
        val effectIdents : jet.MutableSet<String> = TreeSet<String>()
        val humanEffects = ArrayList<Element>()
        val humanInterventions = ArrayList<Element>()
        fun effectIdent(suggestedIdent: String): String{
            // find a unique effect identifier
            var ident = suggestedIdent
            var app = 0
            while (effectIdents.contains(ident)){
                ident = suggestedIdent + app
                app += 1
            }
            return ident
        }
        fun newEffect(ident: String): Element{
            val effect = scenarioDocument.createElement("effect")!!
            if (effectIdents.contains(ident))
                // actually, SetupException isn't really the right choice, but it's not a DocumentException either...
                throw SetupException("effect id is not unique")
            effectIdents.add(ident)
            effect.setAttribute("id",ident)
            humanEffects.add(effect)
            return effect
        }
        fun newIntervention(effects: Array<String>): Element{
            val intervention = scenarioDocument.createElement("intervention")!!
            for (effect in effects){
                val interventionEffect = scenarioDocument.createElement("effect")!!
                interventionEffect.setAttribute("id",effect)
                intervention.appendChild(interventionEffect)
            }
            humanInterventions.add(intervention)
            return intervention
        }
        
        val MDA = getChildElementOpt(interventions, "MDA")
        if (MDA != null){
            val ident = effectIdent("MDA")
            val effect = newEffect(ident)
            val intervention = newIntervention(array(ident))
            
            val timed = getChildElementOpt(MDA,"timed")
            if (timed != null)
                intervention.appendChild(timed) // massList to massListWithCum
            
            effect.appendChild(MDA) // "timed" child removed
        }
        
        val vaccine = getChildElementOpt(interventions, "vaccine")
        if (vaccine != null){
            val ident = effectIdent("vaccine")
            val effect = newEffect(ident)
            val intervention = newIntervention(array(ident))
            val cts = getChildElementOpt(vaccine, "continuous")
            if (cts != null) intervention.appendChild(cts)
            val timed = getChildElementOpt(vaccine, "timed")
            if (timed != null){
                for (deploy in getChildElements(timed, "deploy")){
                    if (deploy.getAttributeNode("cumulativeWithMaxAge") != null)
                        throw DocumentException("can't handle cumulative deployment translations yet")
                }
                intervention.appendChild(timed)
            }
            effect.appendChild(vaccine)
        }
        
        if (humanEffects.size > 0){
            val human = scenarioDocument.createElement("human")!!
            interventions.appendChild(human)
            for (effect in humanEffects)
                human.appendChild(effect)
            for (intervention in humanInterventions)
                human.appendChild(intervention)
        }
    }
}
