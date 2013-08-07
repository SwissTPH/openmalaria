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
import java.util.TreeMap

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
    /** Translate to schema 32.
     * 
     * Many differences to intervention description.
     * 
     * Warning: you may get some validation errors of the type
     * 
     * - Error: org.xml.sax.SAXParseException; cvc-complex-type.2.1: Element
     * - 'cohort' must have no character or element information item
     * - [children], because the type's content type is empty.
     * 
     * This appears to be a bug in the validator, because xmllint reports no
     * problems and the scenarios look valid. */
    public fun translate31To32(){
        var changeIRSReportingToGVI = false
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
        fun newEffect(ident: String, desc: String?): Element{
            val effect = scenarioDocument.createElement("effect")!!
            if (effectIdents.contains(ident))
                // actually, SetupException isn't really the right choice, but it's not a DocumentException either...
                throw SetupException("effect id is not unique")
            effectIdents.add(ident)
            effect.setAttribute("id",ident)
            humanEffects.add(effect)
            if (desc != null)
                effect.setAttribute("name",desc)
            return effect
        }
        
        // move any deployment descriptions from elt to intervention, applying necessary transformations
        fun processDeployments(effects: List<String>, elt: Element, cumCovEffectIdent: String){
            val intervention = scenarioDocument.createElement("intervention")!!
            for (effect in effects){
                val interventionEffect = scenarioDocument.createElement("effect")!!
                interventionEffect.setAttribute("id",effect)
                intervention.appendChild(interventionEffect)
            }
            
            val cts = getChildElementOpt(elt, "continuous")
            if (cts != null) intervention.appendChild(cts)

            val timed = getChildElementOpt(elt, "timed")
            if (timed != null){
                // Note: type change from massList or massCumList to massListWithCum
                // Map of "cumulative max age" to "timed element"
                val cumTimedLists = TreeMap<Double,Element>()
                var hasNonCum = false
                for (deploy in getChildElements(timed, "deploy")){
                    val maxAgeNode = deploy.getAttributeNode("cumulativeWithMaxAge")
                    if (maxAgeNode != null){
                        val age = java.lang.Double.parseDouble( maxAgeNode.getValue() )
                        fun makeTimedList(): Element {
                            val tL = scenarioDocument.createElement("timed")!!
                            val cumCov = scenarioDocument.createElement("cumulativeCoverage")!!
                            cumCov.setAttribute("effect",cumCovEffectIdent) // effect is uniquely used in this case, so translation is exact
                            cumCov.setAttribute("maxAgeYears",java.lang.Double.toString(age))
                            tL.appendChild(cumCov)
                            cumTimedLists.put(age, tL)
                            return tL
                        }
                        val timedList : Element = cumTimedLists.get( age ) ?: makeTimedList()
                        deploy.removeAttribute("cumulativeWithMaxAge")
                        timedList.appendChild(deploy)
                    }else{
                        hasNonCum = true
                    }
                }
                if (hasNonCum)
                    intervention.appendChild(timed)
                else
                    elt.removeChild(timed)
                for (item in cumTimedLists){
                    intervention.appendChild(item.component2())
                }
            }

            val nameAttr = elt.getAttributeNode("name")
            if (nameAttr != null){
                // move name attribute if it exists
                intervention.setAttribute("name",nameAttr.getValue())
                elt.removeAttribute("name")
            }
            
            if (getChildElements(intervention, "continuous").size +
                getChildElements(intervention, "timed").size > 0)
            {
                humanInterventions.add(intervention)
            }
        }
        
        fun updateElt(srcName: String, trgName: String, stripDescElt: Boolean){
            val elt = getChildElementOpt(interventions, srcName)
            if (elt != null){
                val ident = effectIdent(srcName)
                val effect = newEffect(ident, elt.getAttributeNode("name")?.getValue())

                processDeployments(listOf(ident), elt, ident)
                
                if (stripDescElt){
                    val desc = getChildElement(elt,"description")
                    val renamed = scenarioDocument.renameNode(desc,"",trgName)!!
                    effect.appendChild(renamed)
                    interventions.removeChild(elt)  // now defunct
                }else{
                    val renamed = scenarioDocument.renameNode(elt,"",trgName)!!
                    effect.appendChild(renamed) // after removal of "continuous" and "timed" child elements
                }
            }
        }
        fun updateVaccineElt(){
            val elt = getChildElementOpt(interventions, "vaccine")
            if (elt != null){
                val idents = ArrayList<String>()
                for (vacc in getChildElements(elt, "description")){
                    val vaccineType = vacc.getAttribute("vaccineType")!!
                    vacc.removeAttribute("vaccineType")
                    val ident = effectIdent(vaccineType)
                    idents.add(ident)
                    val effect = newEffect(ident, null)
                    val renamed = scenarioDocument.renameNode(vacc,"",vaccineType)!!
                    effect.appendChild(renamed)
                }
                if (idents.size == 0){
                    // corner case: vaccine element with no children was allowed, but is useless
                    interventions.removeChild(elt)
                    return
                }

                // Use any identifier for cumCov since all effects are deployed simultaneously
                processDeployments(idents, elt, idents[0])
                
                interventions.removeChild(elt)  // now defunct
            }
        }
        fun updateIRS(){
            var name = "IRS"
            val elt = getChildElementOpt(interventions, name)
            if (elt != null){
                val desc = getChildElementOpt(elt,"description")
                val renamed = if (desc != null){
                    changeIRSReportingToGVI = true
                    name = "GVI"
                    scenarioDocument.renameNode(desc,"",name)!!
                }else{
                    val desc2 = getChildElement(elt,"description_v2")
                    scenarioDocument.renameNode(desc2,"",name)!!
                }
                val ident = effectIdent(name)
                val effect = newEffect(ident, elt.getAttributeNode("name")?.getValue())
                effect.appendChild(renamed)
                
                processDeployments(listOf(ident), elt, ident)
                interventions.removeChild(elt)  // now defunct
            }
        }
        updateElt("MDA", "MDA", false)
        updateVaccineElt()
        updateElt("IPT", "IPT", true) 
        updateElt("ITN", "ITN", true)
        updateIRS()
        updateElt("vectorDeterrent", "GVI", false)
        updateElt("cohort", "cohort", false)
        
        if (humanEffects.size > 0){
            val human = scenarioDocument.createElement("human")!!
            interventions.appendChild(human)
            for (effect in humanEffects)
                human.appendChild(effect)
            for (intervention in humanInterventions)
                human.appendChild(intervention)
        }
        
        if (changeIRSReportingToGVI){
            val mon = getChildElement(scenarioElement, "monitoring")
            val survOpts = getChildElement(mon, "SurveyOptions")
            for (opt in getChildElements(survOpts, "option")){
                if (opt.getAttribute("name").equals("nMassIRS"))
                    // replace the name (IRS won't be used so don't need both opts)
                    opt.setAttribute("name","nMassGVI")
            }
        }
    }
}
