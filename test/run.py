#!/usr/bin/python
# -*- coding: utf-8 -*-
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

sys.path[0]="@CMAKE_CURRENT_SOURCE_DIR@"
import compareOutputsFloat as compareOuts

# replaced by CMake; run the version it puts in the build/test/ dir.
testSrcDir="@CMAKE_CURRENT_SOURCE_DIR@"
testBuildDir="@CMAKE_CURRENT_BINARY_DIR@"
if not os.path.isdir(testSrcDir) or not os.path.isdir(testBuildDir):
  print "Don't run this script directly; configure CMake then use the version in the CMake build dir."
  sys.exit(-1)
if not os.path.isfile (os.path.join(testSrcDir,"@OM_BOXTEST_SCHEMA_NAME@")):
  print "File not found (wrong CMake var?): ",os.path.join(testSrcDir,"@OM_BOXTEST_SCHEMA_NAME@")
  sys.exit(-1)

# executable
def findFile (*names):
  execs=set()
  for name in names:
    path=os.path.join(testBuildDir,name)
    if os.path.isfile (path):
      execs.add (path)
  
  if not execs:
    print "Unable to find: openMalaria[.exe]; please compile it."
    sys.exit(-1)
  
  newest=None
  for path in execs:
    if newest is None or os.path.getmtime(path) > os.path.getmtime(newest):
       newest=path
  return newest

openMalariaExec=findFile (*["../openMalaria", "../openMalaria.exe", "../debug/openMalaria.exe", "../release/openMalaria.exe"])

def linkOrCopy (src, dest):
  if hasattr(os, 'symlink'):
    os.symlink(os.path.abspath(src), dest)
  else:
    shutil.copy2(src, dest)

# Run, with file "scenario"+name+".xml"
def runScenario(options,omOptions,name):
  scenarioSrc=os.path.join(testSrcDir,"scenario%s.xml" % name)
  if options.xmlValidate:
    return subprocess.call (["xmllint","--noout","--schema","@OM_BOXTEST_SCHEMA_NAME@",scenarioSrc],
			 cwd=testBuildDir)
  
  cmd=options.wrapArgs+[openMalariaExec,"--scenario",scenarioSrc]+omOptions
  
  if not options.run:
    print "\033[1;32m",cmd,"\033[0;00m"
    return 0
  
  # Run from a temporary directory, so checkpoint files won't conflict
  simDir = tempfile.mkdtemp(prefix=name+'-', dir=testBuildDir)
  outFile=os.path.join(simDir,"output.txt")
  checkFile=os.path.join(simDir,"checkpoint")
  
  # Link or copy required files.
  densities_csv=os.path.join(simDir,"densities.csv")
  scenario_xsd=os.path.join(simDir,"@OM_BOXTEST_SCHEMA_NAME@")
  linkOrCopy (os.path.join(testSrcDir,"densities.csv"), densities_csv)
  linkOrCopy (os.path.join(testSrcDir,"@OM_BOXTEST_SCHEMA_NAME@"), scenario_xsd)
  
  if options.logging:
    print time.strftime("\033[0;33m%a, %d %b %Y %H:%M:%S")
  
  lastTime=time.time()
  # While no output and cmd exits successfully:
  while (not os.path.isfile(outFile)):
    if options.logging:
      print "\033[1;32m",cmd,"\033[0;00m"
    ret=subprocess.call (cmd, shell=False, cwd=simDir)
    if ret != 0:
      print "Non-zero exit status: {0}".format(ret)
      break
    
    # if the checkpoint file hasn't been updated, stop
    if not os.path.isfile(checkFile):
      break
    checkTime=os.path.getmtime(checkFile)
    if not checkTime > lastTime:
      break
    lastTime=checkTime
  
  if options.cleanup:
    os.remove(densities_csv)
    os.remove(scenario_xsd)
    for f in (glob.glob(os.path.join(simDir,"checkpoint*")) + glob.glob(os.path.join(simDir,"seed?"))):
      if os.path.isfile(f):
	os.remove(f)
  
  ret = 1
  if os.path.isfile(outFile):
    print "\033[1;34m",
    ret = compareOuts.main (*["",os.path.join(testSrcDir,"original%s.txt"%name), outFile, 1])
    if ret == 0 and options.cleanup:
      os.remove(outFile)
  else:
    stderrFile=os.path.join(simDir,"stderr.txt")
    if os.path.isfile (stderrFile):
      print "\033[0;31mNo results output; error messages:"
      se = file.open(stderrFile)
      se.read()
      se.close()
    else:
      print "\033[0;31mNo results output; error messages:"
  
  try:
    os.rmdir(simDir)
  except OSError:
    print "Directory %s not empty, so not deleted!" % simDir
  
  print "\033[0;00m"
  return ret

def setWrapArgs(option, opt_str, value, parser, *args, **kwargs):
  parser.values.wrapArgs = args[0]

# Test for options
def evalOptions (args):
  parser = OptionParser(usage="Usage: %prog [options] [-- openMalaria options] [scenarios]",
			description="""Scenarios to be run must be of the form scenarioXX.xml; if any are passed on the command line, XX is substituted for each given; if not then all files of the form scenario*.xml are run as test scenarios.
			You can pass options to openMalaria by first specifying -- (to end options passed from the script); for example: %prog 5 -- --print-model""")
  
  parser.add_option("-q","--quiet",
		    action="store_false", dest="logging", default=True,
		    help="Turn off console output from this script")
  parser.add_option("-n","--dry-run", action="store_false", dest="run", default=True,
		    help="Don't actually run openMalaria, just output the commandline.")
  parser.add_option("-c","--dont-cleanup", action="store_false", dest="cleanup", default=True,
		    help="Don't clean up expected files from the temparary dir (checkpoint files, densities.csv and @OM_BOXTEST_SCHEMA_NAME@)")
  parser.add_option("--valid","--validate",
		    action="store_true", dest="xmlValidate", default=False,
		    help="Validate the XML file(s) using xmllint and the latest schema.")
  parser.add_option("-g","--gdb", action="callback", callback=setWrapArgs,
		    callback_args=(["gdb","--args"],),
		    help="Run openMalaria through gdb.")
  parser.add_option("--valgrind", action="callback", callback=setWrapArgs,
		    callback_args=(["valgrind","--gen-suppressions=yes","leak-check=full"],),
		    help="Run openMalaria through valgrind.")
  parser.add_option("--valgrind-track-origins", action="callback", callback=setWrapArgs,
		    callback_args=(["valgrind","--gen-suppressions=yes","leak-check=full","--track-origins=yes"],),
		    help="As --valgrind, but pass --track-origins=yes option (1/2 performance).")
  (options, others) = parser.parse_args(args=args)
  
  options.ensure_value("wrapArgs", [])
  
  toRun=set()
  omOptions=[]
  for arg in others:
    if (arg[0] == '-'):
      omOptions = omOptions + [arg]
    else:
      toRun.add (arg)
  
  return options,omOptions,toRun


def main(args):
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

if __name__ == "__main__":
  sys.exit(main(sys.argv))
