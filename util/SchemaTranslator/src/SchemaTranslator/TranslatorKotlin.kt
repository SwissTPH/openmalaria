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
import org.w3c.dom.Attr

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
enum class HSTreatmentOption {
    NONE; SIMPLE; LEGACY; BLANK
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
    
    var hsTreatmentTranslation = HSTreatmentOption.NONE;
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
            if (schemaVersion < 32){
                scenarioElement.setAttribute("xsi:noNamespaceSchemaLocation", schemaFileName)
            }else{
                // From version 32, we use an explicit namespace:
                scenarioElement.setAttribute("xsi:schemaLocation",
                    "http://openmalaria.org/schema/scenario_" + Integer.toString(schemaVersion)
                    + " " + schemaFileName)
            }

            val translateMeth : String = "translate" + (schemaVersion - 1) + "To" + schemaVersion
            val method : Method? = cls.getMethod(translateMeth)
            if (method == null)
                throw Exception("Method " + translateMeth + " not found")
            try{
                val result = method.invoke(this) // invoke, may throw
                if (result != null){    // null result implies success
                    when (result){
                        // for backwards compatibility:
                        is Boolean -> {
                            if (!result) throw DocumentException("Translation failed (no message)")
                        }
                        else -> throw Exception("Unexpected result value while updating to version ${schemaVersion}:  \"${result}\"")
                    }
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

// For some reason I get type errors if this is declared within the function it's used in
class TimedKey(cohort:Boolean, age:Double?) : Comparable<TimedKey>{
    val cohort = cohort
    val age = age
    override fun compareTo(other: TimedKey): Int =
        if (cohort != other.cohort){
            if (cohort == false) -1 else 1
        }else if (age == null){
            if (other.age == null) 0 else -1
        }else if (other.age == null){
            1
        }else{
            java.lang.Double.compare(age, other.age)
        }
}

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
        var cohortElt : Element? = null
        val interventions = getChildElementOpt(scenarioElement, "interventions")
        if (interventions != null){
            val componentIdents : kotlin.MutableSet<String> = TreeSet<String>()
            val humanComponents = ArrayList<Element>()
            val humanDeployments = ArrayList<Element>()
            fun componentIdent(suggestedIdent: String): String{
                // find a unique component identifier
                var ident = suggestedIdent
                var app = 0
                while (componentIdents.contains(ident)){
                    ident = suggestedIdent + app
                    app += 1
                }
                return ident
            }
            fun newComponent(ident: String, desc: String?): Element{
                val component = scenarioDocument.createElement("component")!!
                if (componentIdents.contains(ident))
                    // actually, SetupException isn't really the right choice, but it's not a DocumentException either...
                    throw SetupException("component id is not unique")
                componentIdents.add(ident)
                component.setAttribute("id",ident)
                humanComponents.add(component)
                if (desc != null)
                    component.setAttribute("name",desc)
                return component
            }
            
            // move any deployment descriptions from elt to deployment, applying necessary transformations
            fun processDeployments(components: List<String>, elt: Element, cumCovComponentId: String,
                addVaccLimits: Boolean){
                val deployment = scenarioDocument.createElement("deployment")!!
                for (component in components){
                    val deploymentComponent = scenarioDocument.createElement("component")!!
                    deploymentComponent.setAttribute("id",component)
                    deployment.appendChild(deploymentComponent)
                }
                
                fun extractBool(attr: Attr?, deploy: Element): Boolean {
                    if (attr == null) return false
                    else{
                        val result = java.lang.Boolean.parseBoolean(attr.getValue())
                        deploy.removeAttributeNode(attr)
                        return result
                    }
                }

                val cts = getChildElementOpt(elt, "continuous")
                if (cts != null) {
                    var ctsNonCohort: Element? = null
                    var ctsCohort: Element? = null
                    
                    var n : Int = 0
                    for (deploy in getChildElements(cts, "deploy")){
                        if (addVaccLimits){
                            deploy.setAttribute("vaccMinPrevDoses", Integer.toString(n))
                            deploy.setAttribute("vaccMaxCumDoses", Integer.toString(n+1))
                        }
                        
                        val cohort : Boolean = extractBool(deploy.getAttributeNode("cohort"), deploy)
                        val list: Element = if (cohort){
                            if (ctsCohort == null) {
                                ctsCohort = scenarioDocument.createElement("continuous")!!
                                val rTC = scenarioDocument.createElement("restrictToCohort")!!
                                rTC.setAttribute("id","cohort") // component id _will_ be "cohort"
                                ctsCohort!!.appendChild(rTC)
                            }
                            ctsCohort!!
                        }else{
                            if (ctsNonCohort == null) {
                                ctsNonCohort = scenarioDocument.createElement("continuous")!!
                            }
                            ctsNonCohort!!
                        }
                        list.appendChild(deploy)
                        n += 1
                    }
                    
                    elt.removeChild(cts)
                    if (ctsNonCohort != null) deployment.appendChild(ctsNonCohort!!)
                    if (ctsCohort != null) deployment.appendChild(ctsCohort!!)
                }

                val timed = getChildElementOpt(elt, "timed")
                if (timed != null){
                    // Note: type change from massList or massCumList to massListWithCum
                    // Map of "cumulative max age" to "timed element"
                    val cumTimedLists = TreeMap<TimedKey,Element>()
                    for (deploy in getChildElements(timed, "deploy")){
                        fun extractDoubleOpt(attr: Attr?): Double? {
                            if (attr == null) return null
                            else{
                                val result = java.lang.Double.parseDouble(attr.getValue())
                                deploy.removeAttributeNode(attr)
                                return result
                            }
                        }

                        val cohort : Boolean = extractBool(deploy.getAttributeNode("cohort"), deploy)
                        val cumAge : Double? = extractDoubleOpt(deploy.getAttributeNode("cumulativeWithMaxAge"))
                        val key = TimedKey(cohort, cumAge)

                        fun makeTimedList(): Element {
                            val tL = scenarioDocument.createElement("timed")!!
                            if (cohort){
                                val rTC = scenarioDocument.createElement("restrictToCohort")!!
                                rTC.setAttribute("id","cohort") // component id _will_ be "cohort"
                                tL.appendChild(rTC)
                            }
                            if (cumAge != null){
                                val cumCov = scenarioDocument.createElement("cumulativeCoverage")!!
                                cumCov.setAttribute("component",cumCovComponentId) // component is uniquely used in this case, so translation is exact
                                cumCov.setAttribute("maxAgeYears",java.lang.Double.toString(cumAge))
                                tL.appendChild(cumCov)
                            }
                            cumTimedLists.put(key, tL)
                            return tL
                        }
                        val timedList : Element = cumTimedLists.get( key ) ?: makeTimedList()
                        deploy.removeAttribute("cumulativeWithMaxAge")
                        timedList.appendChild(deploy)
                    }
                    elt.removeChild(timed)
                    for (item in cumTimedLists){
                        deployment.appendChild(item.component2())
                    }
                }

                val nameAttr = elt.getAttributeNode("name")
                if (nameAttr != null){
                    // move name attribute if it exists
                    deployment.setAttribute("name",nameAttr.getValue())
                    elt.removeAttribute("name")
                }
                
                if (getChildElements(deployment, "continuous").size +
                    getChildElements(deployment, "timed").size > 0)
                {
                    humanDeployments.add(deployment)
                }
            }
            
            fun updateElt(srcName: String, trgName: String, stripDescElt: Boolean): Element?{
                val elt = getChildElementOpt(interventions, srcName)
                if (elt != null){
                    val ident = componentIdent(srcName)
                    val component = newComponent(ident, elt.getAttributeNode("name")?.getValue())

                    processDeployments(listOf(ident), elt, ident, false)
                    
                    if (stripDescElt){
                        val desc = getChildElement(elt,"description")
                        val renamed = scenarioDocument.renameNode(desc,"",trgName)!!
                        component.appendChild(renamed)
                        interventions.removeChild(elt)  // now defunct
                    }else{
                        val renamed = scenarioDocument.renameNode(elt,"",trgName)!!
                        component.appendChild(renamed) // after removal of "continuous" and "timed" child elements
                    }
                    return component
                }
                return null
            }
            fun updateMDA(){
                val name = "MDA"
                val elt = getChildElementOpt(interventions, name)
                if (elt != null){
                    val ident = componentIdent(name)
                    val component = newComponent(ident, elt.getAttributeNode("name")?.getValue())

                    processDeployments(listOf(ident), elt, ident, false)
                    
                    val desc1d = getChildElementOpt(elt, "description")
                    if (desc1d != null ){
                        val renamed = scenarioDocument.renameNode(desc1d, "", "MDA1D")!!
                        component.appendChild(renamed)
                        interventions.removeChild(elt)
                    }else{
                        // no 1-day-TS description; add the new 5-day-TS drug description
                        val effects = scenarioDocument.createElement("effects")!!
                        elt.appendChild(effects)
                        val option = scenarioDocument.createElement("option")!!
                        effects.appendChild(option)
                        option.setAttribute("name","clear blood-stage infections")
                        option.setAttribute("pSelection","1")
                        val clearInfections = scenarioDocument.createElement("clearInfections")!!
                        option.appendChild(clearInfections)
                        clearInfections.setAttribute("timesteps","1")
                        clearInfections.setAttribute("stage","blood")
                        component.appendChild(elt) // after removal of "continuous" and "timed" child elements
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
                        val ident = componentIdent(vaccineType)
                        idents.add(ident)
                        val component = newComponent(ident, null)
                        val renamed = scenarioDocument.renameNode(vacc,"",vaccineType)!!
                        component.appendChild(renamed)
                    }
                    if (idents.size == 0){
                        // corner case: vaccine element with no children was allowed, but is useless
                        interventions.removeChild(elt)
                        return
                    }

                    // Use any identifier for cumCov since all components are deployed simultaneously
                    processDeployments(idents, elt, idents[0], true)
                    
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
                    val ident = componentIdent(name)
                    val component = newComponent(ident, elt.getAttributeNode("name")?.getValue())
                    component.appendChild(renamed)
                    
                    processDeployments(listOf(ident), elt, ident, false)
                    interventions.removeChild(elt)  // now defunct
                }
            }
            
            updateMDA()
            updateVaccineElt()
            updateElt("IPT", "IPT", true) 
            updateElt("ITN", "ITN", true)
            updateIRS()
            updateElt("vectorDeterrent", "GVI", false)
            cohortElt = updateElt("cohort", "cohort", false)
            updateElt("immuneSuppression", "clearImmunity", false)
            
            if (humanComponents.size > 0){
                val human = scenarioDocument.createElement("human")!!
                interventions.appendChild(human)
                for (component in humanComponents)
                    human.appendChild(component)
                for (deployment in humanDeployments)
                    human.appendChild(deployment)
            }
        }

        fun updateHSIO(immOut: Element){
            val drugReg = getChildElement(immOut, "drugRegimen")
            val drugs = TreeSet<String>()
            drugs.add(drugReg.getAttribute("firstLine")!!)
            drugs.add(drugReg.getAttribute("secondLine")!!)
            drugs.add(drugReg.getAttribute("inpatient")!!)
            val treatActions = scenarioDocument.createElement("treatmentActions")!!
            for (drug in drugs){
                val x = scenarioDocument.createElement(drug)!!
                treatActions.appendChild(x)     //TODO: order
                if (options.hsTreatmentTranslation == HSTreatmentOption.NONE){
                    throw DocumentException("please re-run, setting the --hsTreatmentTranslation option")
                }else if (options.hsTreatmentTranslation == HSTreatmentOption.BLANK){
                    x.appendChild(scenarioDocument.createComment("<clearInfections .../>")!!)
                }else if (options.hsTreatmentTranslation == HSTreatmentOption.LEGACY){
                    x.setAttribute("name","legacy (emulate pre-32 treatment)")
                    val clear = scenarioDocument.createElement("clearInfections")!!
                    clear.setAttribute("timesteps","-1")
                    clear.setAttribute("stage","both")
                    x.appendChild(clear)
                }else if (options.hsTreatmentTranslation == HSTreatmentOption.SIMPLE){
                    x.setAttribute("name","clear blood-stage infections")
                    val clear = scenarioDocument.createElement("clearInfections")!!
                    clear.setAttribute("timesteps","1")
                    clear.setAttribute("stage","blood")
                    x.appendChild(clear)
                }
            }
            immOut.insertBefore(treatActions,getChildElement(immOut,"pSeekOfficialCareUncomplicated1"))
        }
        
        val monitoring = getChildElementOpt(scenarioElement, "monitoring")
        if (monitoring != null){
            if (changeIRSReportingToGVI){
                val survOpts = getChildElement(monitoring, "SurveyOptions")
                for (opt in getChildElements(survOpts, "option")){
                    if (opt.getAttribute("name").equals("nMassIRS"))
                        // replace the name (IRS won't be used so don't need both opts)
                        opt.setAttribute("name","nMassGVI")
                }
            }
            
            if (cohortElt != null){
                fun copyAttr(attrName: String){
                    val attr = monitoring.getAttributeNode(attrName)
                    if (attr != null){
                        cohortElt!!.setAttribute(attrName, attr.getValue())
                    }
                }
                copyAttr("firstBoutOnly")
                copyAttr("firstTreatmentOnly")
                copyAttr("firstInfectionOnly")
            }
            // these we need to do even if cohortElt == null:
            monitoring.removeAttribute("firstBoutOnly")
            monitoring.removeAttribute("firstTreatmentOnly")
            monitoring.removeAttribute("firstInfectionOnly")

            val survOpts = getChildElement(monitoring, "SurveyOptions")
            for (opt in getChildElements(survOpts, "option")){
                if (opt.getAttribute("name") == "nMassVA")
                    opt.setAttribute("name", "nMassGVI")
            }
        }
        
        val healthSystem = getChildElementOpt(scenarioElement, "healthSystem")
        if (healthSystem != null){
            val immOut = getChildElementOpt(healthSystem, "ImmediateOutcomes")
            if (immOut != null) updateHSIO( immOut )
        }
        if (interventions != null){
            val changeHS = getChildElementOpt(interventions, "changeHS")
            if (changeHS != null){
                for (deploys in getChildElements(changeHS, "timedDeployment")){
                    val immOut = getChildElementOpt(deploys, "ImmediateOutcomes")
                    if (immOut != null) updateHSIO( immOut )
                }
            }
        }
        
        val model = getChildElement(scenarioElement, "model")
        val modelOpts = getChildElement(model, "ModelOptions")
        val opt = scenarioDocument.createElement("option")!!
        opt.setAttribute("name", "INDIRECT_MORTALITY_FIX")
        // This option could affect parameterisation, best to leave the bug fix off by default
        opt.setAttribute("value", "false")
        modelOpts.appendChild(opt)
        
        // From version 32, we use an explicit namespace; delete the old here,
        // the new is set by translateAndValidate()
        scenarioElement.removeAttribute("xsi:noNamespaceSchemaLocation")
        scenarioDocument.renameNode(scenarioElement, 
            "http://openmalaria.org/schema/scenario_32", "om:scenario")
        scenarioElement = scenarioDocument.getDocumentElement()!!
    }
    
    /** Translate to schema 33.
     *
     * TODO: also "screen"
     * 
     * MDA descriptions changed as did 1-day TS health system configurations
     * (the latter will not be automatically updated). Many other
     * backwards-compatible changes occurred. */
    public fun translate32To33(){
        throw DocumentException("This function is incomplete. Suggest using translateXML.py instead.")
        if (getChildElementOpt(getChildElement(scenarioElement, "healthSystem"), "EventScheduler") != null ){
            throw DocumentException("Refusing to translate EventScheduler element to schema 33")
        }
        val interventions = getChildElementOpt(scenarioElement, "interventions")
        if( interventions != null ){
            val human = getChildElementOpt(interventions, "human")
            if( human != null ){
                for (comp in getChildElements(human, "component")){
                    // Component MDA: replace with treatSimple
                    val mda = getChildElementOpt(comp, "MDA")
                    if( mda != null ){
                        val eff = getChildElement(mda, "effects")
                        val opts = getChildElements(eff, "option")
                        if( opts.size == 1 ){
                            val opt = opts[0]
                            if( java.lang.Double.parseDouble(opt.getAttribute("pSelection")) != 1.0 ){
                                throw DocumentException("Invalid document: MDA effects has single option, pSelection != 1")
                            }
                            if( getChildElementOpt(opt, "deploy") != null ){
                                throw DocumentException("Lazily refusing to translate triggered deployment in MDA")
                            }
                            var durLiver = 0
                            var durBlood = 0  // durations, in steps
                            for (clearInf in getChildElements(opt, "clearInfections")){
                                val ts = Integer.parseInt(clearInf.getAttribute("timesteps"))
                                when (clearInf.getAttribute("stage")){
                                    "liver" -> durLiver = ts
                                    "blood" -> durBlood = ts
                                    "both" -> {
                                        durLiver = ts
                                        durBlood = ts
                                    }
                                    else -> throw DocumentException("Invalid document: clearInfections' stage attribute should be liver, blood or both")
                                }
                            }
                            val treatSimple = scenarioDocument.createElement("treatSimple")!!
                            treatSimple.setAttribute("durationLiver", "${durLiver}t")
                            treatSimple.setAttribute("durartionBlood", "${durBlood}t")
                            comp.removeChild(opt)
                            comp.appendChild(treatSimple)
                        }else{
                            throw DocumentException("Lazily refusing to translate MDA with number options > 1")
                        }
                    }else if( getChildElementOpt(comp, "MDA1D") != null ){
                        throw DocumentException("Refusing to translate MDA1D element to schema 33")
                    }
                }
            }
        }
    }
}
