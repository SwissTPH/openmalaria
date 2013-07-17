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
import java.io.FileReader
import java.io.FileOutputStream
import javax.xml.transform.stream.StreamResult
import javax.xml.transform.TransformerFactory
import javax.xml.transform.Result
import javax.xml.transform.OutputKeys
import javax.xml.transform.dom.DOMSource
import javax.xml.validation.SchemaFactory
import javax.xml.validation.Schema
import javax.xml.validation.Validator
import javax.xml.transform.Source
import org.xml.sax.SAXParseException
import java.sql.Connection
import java.sql.DriverManager
import java.sql.Statement
import java.sql.ResultSet
import java.io.StringReader
import java.io.StringWriter
import java.lang.reflect.Method


/** Entrypoint for the SchemaTranslator. This depends on both TranslatorKotlin and TranslatorJava. */

fun main(args: Array<String>): Unit {
    val options = Options()
    try{
        processArgs (args, options)
    }catch(e: SetupException){
        System.err.println("""Error: ${e.getMessage()}

Usage: schemaTranslator [options]:
  --required-version VERSION	The version number to update the document(s) to.
				Default is latest version: ${options.latestVersion}.
  --latest-schema		Use schema scenario.xsd instead of scenario_XX.xsd
  --current-schema		Use schema scenario_current.xsd instead of scenario_XX.xsd
  --no-validation		Don't validate the result
  --no-translation		Don't write out the translated result (but still
				translate internally for validation)
  --update-db			Update DB entries instead of files
  --maxDensCorrection BOOL	Update 12->13 requires this sometimes: set true to
				include bug fix, false to explicitly exclude it.
  --iptiSpOptionWithoutInterventions
				For scenarios with iptiDescription but without
				interventions, assume usage of the IPTI model
				was (t) intended or (f) a mistake.
  --iptiReportOnlyAtRisk BOOL	Previously the IPTI_SP_MODEL option implied this
				behaviour although not necessarily intended; now
				this behaviour is controlled by the separate option
				REPORT_ONLY_AT_RISK. Specifying true here causes
				option REPORT_ONLY_AT_RISK to be added to scenarios
				already using IPTI_SP_MODEL.
  --ITN-description ARG		Update 29: new parameterisation is needed.
				Use ARG "replace" to replace old parameters with a new default
				parameterisation, or "manual" if you will update by hand
				(in the second case, you should use --no-validation).
  --schema-folder		The schema folder, by default ../../schema
  --input-folder		The input folder, by default ./scenarios
  --output-folder		The output folder, by default ./translatedScenarios
""")
        System.exit(1)
    }

    try{
        if (options.doDBUpdate){
            updateDB(options)
        }else{
            System.out.println("Put XMLs to be translated into the \"${options.inputFolder}\" directory")
            visitAllFiles(options)
        }
    }catch (e : Exception) {
        e.printStackTrace()
        System.exit(1)
    }
}

fun processArgs(args: Array<String>, options: Options) {
    var i : Int = 0
    while (i < args.size){
        var arg : String = args[i].replace('_', '-')
        when (arg){
            "--required-version" -> options.targetVersion = Integer.parseInt(args[++i])
            "--oneDayTimesteps" -> {
                options.doODTTranslation = true
                options.doValidation = false
                System.out.println("You have chosen the --oneDayTimesteps option, this option is only intended for the fitting scenarii or scenarii using no intervention. Validation disabled.")
            }
            "--latest-schema" -> {
                if (options.latestSchema != SchemaName.VERSIONED){
                    System.err.println("Warning: overwriting another schema-name option")
                }
                options.latestSchema = SchemaName.NO_SUFFIX
            }
            "--current-schema" -> {
                if (options.latestSchema != SchemaName.VERSIONED){
                    System.err.println("Warning: overwriting another schema-name option")
                }
                options.latestSchema = SchemaName.CURRENT
            }
            "--no-validation" -> options.doValidation = false
            "--no-translation" -> options.doTranslation = false
            "--update-db" -> options.doDBUpdate = true
            "--maxDensCorrection" -> {
                options.maxDensBug = when (args[++i].toLowerCase()){
                    "true" -> BugCorrectionBehaviour.CORRECT
                    "false" -> BugCorrectionBehaviour.DONT_CORRECT
                    else -> throw SetupException("--maxDensCorrection: expected true or false")
                }
            }
            "--iptiSpOptionWithoutInterventions" -> {
                options.iptiSpOption = when (args[++i].toLowerCase()){
                    "true" -> IptiSpBehaviour.ASSUME_INTENDED
                    "false" -> IptiSpBehaviour.ASSUME_UNINTENDED
                    else -> throw SetupException("--iptiSpOptionWithoutInterventions: expected true or false")
                }
            }
            "--iptiReportOnlyAtRisk" -> {
                options.iptiROAR = when (args[++i].toLowerCase()){
                    "true" -> IptiReportOnlyAtRiskBehaviour.ON
                    "false" -> IptiReportOnlyAtRiskBehaviour.OFF
                    else -> throw SetupException("--iptiReportOnlyAtRisk: expected true or false")
                }
            }
            "--schema-folder" -> options.schemaFolder = File(args[++i])
            "--input-folder" -> options.inputFolder = File(args[++i])
            "--output-folder" -> options.outputFolder = File(args[++i])
            "--ITN-description" -> {
                if( options.ITN29Translation != ITN29ParameterTranslation.NONE )
                    throw SetupException("--ITN-description given multiple times")
                options.ITN29Translation = when (args[++i].toLowerCase()){
                    "replace" -> ITN29ParameterTranslation.REPLACE
                    "manual" -> ITN29ParameterTranslation.MANUAL
                    else -> throw SetupException("--ITN-description received unexpected argument: ${args[i]}")
                }
            }
            else -> throw SetupException("unknown option: ${arg}")
        }
        i++
    }

    if (options.targetVersion == 1)
        throw SetupException("Target version 1 is not supported (see comment for translate0To1).")
    checkFolder(options.inputFolder,true)
    checkFolder(options.outputFolder,true)
    checkFolder(options.schemaFolder,false)
}

fun checkFolder(file: File, create: Boolean = false){
    if (!file.exists()){
        if (create){
            try{
                file.mkdir()
            }catch(e: SecurityException){
                throw SetupException("unable to create directory ${file.path}")
            }
        }else{
            throw SetupException("no such directory: ${file.path}")
        }
    }
    if (!file.isDirectory()){
        throw SetupException("${file.path} is not a directory")
    }
}

fun visitAllFiles(options: Options) : Unit {
    var errorCount = 0
    fun visit(input: File, outDir: File, top: Boolean){
        if (input.isDirectory()){
            var subDir = if (top) outDir else File(outDir, input.getName())
            for( child in (input.list() as Array<String>)){
                visit(File(input, child), subDir, false)
            }
        }else{
            if (input.getName().endsWith(".xml")){
                System.out.println("Translating ${input.getAbsolutePath()}")
                try{
                    val translator = TranslatorJava(InputSource(FileReader(input)), options)
                    translator.translateAndValidate()
                    if (options.doTranslation||options.doODTTranslation) {
                        val outFile = File(outDir, input.getName())
                        outFile.createNewFile()
                        val result = StreamResult(FileOutputStream(outFile))
                        translator.writeTo (result)
                    }
                }catch(e: DocumentException){
                    System.err.println("Error updating a document: ${e.getMessage()}")
                    if (errorCount < 2) e.printStackTrace()
                    ++errorCount
                    if (errorCount > 10){
                        throw Exception("Too many errors encountered. Stopping.")
                    }
                }
            }
        }
    }
    visit(options.inputFolder, options.outputFolder, true)
}

fun updateDB(options: Options) {
    Class.forName("com.mysql.jdbc.Driver").newInstance()
    val con: Connection = DriverManager.getConnection("jdbc:mysql://127.0.0.1:3306/DBNAME", "USER", "PASSWD")
    try{
        if (con.isClosed()) {
            System.err.println("Failed to connect to MySQL server")
            return
        }
        System.err.println("Successfully connected to MySQL server using TCP/IP...")

        val stmt : Statement = con.createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_UPDATABLE)!!
        val uprs : ResultSet = stmt.executeQuery("SELECT id,run_id,xml FROM scenarios")
        var errorCount = 0
        while (uprs.next()){
            val xml : String = uprs.getString("xml")!!
            val inputSource : InputSource = InputSource()
            inputSource.setCharacterStream(StringReader(xml))
            try{
                val translator = TranslatorJava(inputSource, options)
                translator.translateAndValidate()
                val result : StreamResult = StreamResult(StringWriter())
                translator.writeTo (result)
                val translatedXML : String = result.getWriter()!!.toString()!!
                uprs.updateString("xml", translatedXML)
                            
                uprs.updateRow()
            }catch(e: DocumentException){
                System.err.println("Error updating a document: ${e.getMessage()}")
                if (errorCount < 2) e.printStackTrace()
                ++errorCount
                if (errorCount > 10){
                    throw Exception("Too many errors encountered. Stopping.")
                }
            }
        }
    }finally{
        con.close()
    }
}
