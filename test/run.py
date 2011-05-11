#!/usr/bin/python
# -*- coding: utf-8 -*-

# This file is part of OpenMalaria.
# 
# Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# 
# OpenMalaria is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

# With no arguments, run all scenario*.xml files.
# With arguments A,...,Z, only run scenarioA.xml, ..., scenarioZ.xml
# Exit status:
#	0 - all tests passed (or no tests)
#	1 - a test failed
#	-1 - unable to run test

import sys
import os
import tempfile
import glob
import time
import subprocess
import shutil
from optparse import OptionParser
import gzip

sys.path[0]="@CMAKE_SOURCE_DIR@/util"
import compareOutput
import compareCtsout
import xml.sax.handler

class RunError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

# replaced by CMake; run the version it puts in the build/test/ dir.
testSrcDir="@CMAKE_CURRENT_SOURCE_DIR@"
testBuildDir="@CMAKE_CURRENT_BINARY_DIR@"
if not os.path.isdir(testSrcDir) or not os.path.isdir(testBuildDir):
    print "Don't run this script directly; configure CMake then use the version in the CMake build dir."
    sys.exit(-1)

# executable
def findFile (*names):
    execs=set()
    for name in names:
        path=os.path.join(testBuildDir,name)
        if os.path.isfile (path):
            execs.add (path)
    
    if not execs:
        raise RunError("Unable to find: openMalaria[.exe]; please compile it.")
    
    newest=None
    for path in execs:
        if newest is None or os.path.getmtime(path) > os.path.getmtime(newest):
            newest=path
            return newest
 
class ElementHandler(xml.sax.handler.ContentHandler):
    def ElementHandler():
        self.schema = None
    def startElement(self, name, attributes):
        if name == "scenario":
            try:
                self.schema = attributes["xsi:noNamespaceSchemaLocation"]
            except KeyError:
                raise RunError("can't find noNamespaceSchemaLocation attribute")

def getSchemaName(scenarioFile):
    parser = xml.sax.make_parser()
    handler = ElementHandler()
    parser.setContentHandler(handler)
    parser.parse(scenarioFile)
    if handler.schema == None:
        raise RunError("no scenario element found??")
    return handler.schema

openMalariaExec=os.path.abspath(findFile (*["../openMalaria", "../Debug/openMalaria", "../Release/openMalaria", "../openMalaria.exe", "../debug/openMalaria.exe", "../release/openMalaria.exe"]))

def linkOrCopy (src, dest):
    if not os.path.isfile(src):
        raise RunError("linkOrCopy: can't find file "+src)
    if hasattr(os, 'symlink'):
        os.symlink(os.path.abspath(src), dest)
    else:
        shutil.copy2(src, dest)

# Run, with file "scenario"+name+".xml"
def runScenario(options,omOptions,name):
    scenarioSrc=os.path.abspath(os.path.join(testSrcDir,"scenario%s.xml" % name))
    if not os.path.isfile(scenarioSrc):
        raise RunError('No such scenario file '+scenarioSrc+'!')
    schemaName=getSchemaName(scenarioSrc)
    scenarioSchema=os.path.abspath(os.path.join(testSrcDir,'../schema',schemaName))
    if not os.path.isfile(scenarioSchema):
        scenarioSchema=os.path.abspath(os.path.join(testBuildDir,'../schema',schemaName))
        if not os.path.isfile(scenarioSchema):
            raise RunError("can't find "+schemaName)
    if options.xmlValidate:
        cmd=["xmllint","--noout","--schema",scenarioSchema,scenarioSrc]
        # alternative: ["xmlstarlet","val","-s",SCHEMA,scenarioSrc]
        if options.logging:
            print "\033[0;32m  "+(" ".join(cmd))+"\033[0;00m"
        return subprocess.call (cmd,cwd=testBuildDir)
    
    cmd=options.wrapArgs+[openMalariaExec,"--resource-path",os.path.abspath(testSrcDir),"--scenario",scenarioSrc]+omOptions
    
    if not options.run:
        print "\033[0;32m  "+(" ".join(cmd))+"\033[0;00m"
        return 0
    
    # Run from a temporary directory, so checkpoint files won't conflict
    simDir = tempfile.mkdtemp(prefix=name+'-', dir=testBuildDir)
    outputFile=os.path.join(simDir,"output.txt")
    outputGzFile=os.path.join(simDir,"output.txt.gz")
    ctsoutFile=os.path.join(simDir,"ctsout.txt")
    ctsoutGzFile=os.path.join(simDir,"ctsout.txt.gz")
    checkFile=os.path.join(simDir,"checkpoint")
    
    # Link or copy required files.
    # The schema file only needs to be copied in BOINC mode, since otherwise the
    # scenario is opened with a path and the schema can be found in the same
    # directory. We copy it anyway.
    scenario_xsd=os.path.join(simDir,schemaName)
    linkOrCopy (scenarioSchema, scenario_xsd)
    
    if options.logging:
        print time.strftime("\033[0;33m%a, %d %b %Y %H:%M:%S")+"\t\033[1;33mscenario%s.xml" % name
    
    startTime=lastTime=time.time()
    # While no output.txt file and cmd exits successfully:
    while (not os.path.isfile(outputFile)):
        if options.logging:
            print "\033[0;32m  "+(" ".join(cmd))+"\033[0;00m"
        ret=subprocess.call (cmd, shell=False, cwd=simDir)
        if ret != 0:
            print "\033[1;31mNon-zero exit status: " + str(ret)
            break
        
        # check for output.txt.gz in place of output.txt and uncompress:
        if (os.path.isfile(outputGzFile)) and (not os.path.isfile(outputFile)):
            f_in = gzip.open(outputGzFile, 'rb')
            f_out = open(outputFile, 'wb')
            f_out.writelines(f_in)
            f_out.close()
            f_in.close()
            os.remove(outputGzFile)
        
        # check for ctsout.txt.gz in place of ctsout.txt and uncompress:
        if (os.path.isfile(ctsoutGzFile)) and (not os.path.isfile(ctsoutFile)):
            f_in = gzip.open(ctsoutGzFile, 'rb')
            f_out = open(ctsoutFile, 'wb')
            f_out.writelines(f_in)
            f_out.close()
            f_in.close()
            os.remove(ctsoutGzFile)
        
        # if the checkpoint file hasn't been updated, stop
        if not os.path.isfile(checkFile):
            break
        checkTime=os.path.getmtime(checkFile)
        if not checkTime > lastTime:
            break
        lastTime=checkTime
    
    if ret == 0 and options.logging:
        print "\033[0;33mDone in " + str(time.time()-startTime) + " seconds"
    
    if options.cleanup:
        os.remove(scenario_xsd)
        for f in (glob.glob(os.path.join(simDir,"checkpoint*")) + glob.glob(os.path.join(simDir,"seed?")) + [os.path.join(simDir,"init_data.xml"),os.path.join(simDir,"boinc_finish_called"),os.path.join(simDir,"scenario.sum")]):
            if os.path.isfile(f):
                os.remove(f)
    
    # Compare outputs:
    if ret == 0:
        # ctsout.txt (this output is optional):
        if os.path.isfile(ctsoutFile):
            origCtsout = os.path.join(testSrcDir,"expected/ctsout%s.txt"%name)
            newCtsout = os.path.join(testBuildDir,"ctsout%s.txt"%name)
            if os.path.isfile(origCtsout):
                ctsret,ctsident = compareCtsout.main (origCtsout, ctsoutFile)
            else:
                ctsret,ctsident = 3,False
                print "\033[1;31mNo original ctsout.txt to compare with."
            if ctsident and options.cleanup:
                os.remove(ctsoutFile)
                if os.path.isfile(newCtsout):
                    os.remove(newCtsout)
            else:
                shutil.copy2(ctsoutFile, newCtsout)
                if options.diff:
                    subprocess.call (["kdiff3",origCtsout,ctsoutFile])
        else:
            ctsret,ctsident = 0,True
        
        # output.txt (this output is required):
        if os.path.isfile(outputFile):
            origOutput = os.path.join(testSrcDir,"expected/output%s.txt"%name)
            newOutput = os.path.join(testBuildDir,"output%s.txt"%name)
            if os.path.isfile(origOutput):
                ret,ident = compareOutput.main (origOutput, outputFile, 0)
            else:
                ret,ident = 3,False
                print "\033[1;31mNo original output.txt to compare with."
            if ident and options.cleanup:
                os.remove(outputFile)
                if os.path.isfile(newOutput):
                    os.remove(newOutput)
            else:
                shutil.copy2(outputFile, newOutput)
                if options.diff:
                    subprocess.call (["kdiff3",origOutput,outputFile])
        else:
            ret,ident = 1,False
            stderrFile=os.path.join(simDir,"stderr.txt")
            if os.path.isfile (stderrFile):
                print "\033[1;31mNo output 'output.txt'; error messages:"
                se = open(stderrFile)
                se.read()
                se.close()
            else:
                print "\033[1;31mNo output 'output.txt'"
        
        ret=max(ret,ctsret)
        ident=ident and ctsident
    
    try:
        os.rmdir(simDir)
    except OSError:
        print "\033[0;31mDirectory %s not empty, so not deleted!" % simDir
    
    print "\033[0;00m"
    return ret

def setWrapArgs(option, opt_str, value, parser, *args, **kwargs):
    parser.values.wrapArgs = args[0]

# Test for options
def evalOptions (args):
    #First separate OpenMalaria args and args for this script
    omArgsBegin = len(args)
    for i in range(0,len(args)-1):
        if args[i] == "--":
            omArgsBegin = i+1
            break
    omOptions=args[omArgsBegin:]
    args = args[:omArgsBegin]
    
    parser = OptionParser(usage="Usage: %prog [options] [-- openMalaria options] [scenarios]",
			description="""Scenarios to be run must be of the form scenarioXX.xml; if any are passed on the command line, XX is substituted for each given; if not then all files of the form scenario*.xml are run as test scenarios.
You can pass options to openMalaria by first specifying -- (to end options passed from the script); for example: %prog 5 -- --print-model""")
    
    parser.add_option("-q","--quiet",
		    action="store_false", dest="logging", default=True,
		    help="Turn off console output from this script")
    parser.add_option("-n","--dry-run", action="store_false", dest="run", default=True,
		    help="Don't actually run openMalaria, just output the commandline.")
    parser.add_option("-c","--dont-cleanup", action="store_false", dest="cleanup", default=True,
		    help="Don't clean up expected files from the temparary dir (checkpoint files, schema, etc.)")
    parser.add_option("-d","--diff", action="store_true", dest="diff", default=False,
            help="Launch a diff program (kdiff3) on the output if validation fails")
    parser.add_option("--valid","--validate",
		    action="store_true", dest="xmlValidate", default=False,
		    help="Validate the XML file(s) using xmllint and the latest schema.")
    parser.add_option("-g","--gdb", action="callback", callback=setWrapArgs,
		    callback_args=(["gdb","--args"],),
		    help="Run openMalaria through gdb.")
    parser.add_option("--valgrind", action="callback", callback=setWrapArgs,
		    callback_args=(["valgrind","--gen-suppressions=yes","--leak-check=full"],),
		    help="Run openMalaria through valgrind (using memcheck tool).")
    parser.add_option("--valgrind-track-origins", action="callback", callback=setWrapArgs,
		    callback_args=(["valgrind","--gen-suppressions=yes","--leak-check=full","--track-origins=yes"],),
		    help="As --valgrind, but pass --track-origins=yes option (1/2 performance).")
    parser.add_option("--callgrind", action="callback", callback=setWrapArgs,
            callback_args=(["valgrind","--tool=callgrind"],),
            help="Run openMalaria through valgrind using callgrind tool.")
    parser.add_option("--cachegrind", action="callback", callback=setWrapArgs,
            callback_args=(["valgrind","--tool=cachegrind"],),
            help="Run openMalaria through valgrind using cachegrind tool.")
    (options, others) = parser.parse_args(args=args)
    
    options.ensure_value("wrapArgs", [])
    
    toRun=set()
    for arg in others:
        toRun.add (arg)
    
    return options,omOptions,toRun


def main(args):
    try:
        (options,omOptions,toRun) = evalOptions (args[1:])
        
        if not toRun:
            for p in glob.iglob(os.path.join(testSrcDir,"scenario*.xml")):
                f = os.path.basename(p)
                n=f[8:-4]
                assert ("scenario%s.xml" % n) == f
                toRun.add(n)
        
        retVal=0
        for name in toRun:
            r=runScenario(options,omOptions,name)
            retVal = r if retVal == 0 else retVal
        
        return retVal
    except RunError,e:
        print str(e)
        return -1

if __name__ == "__main__":
    sys.exit(main(sys.argv))
