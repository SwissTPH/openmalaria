# -*- coding: utf-8 -*-
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

# This is a minimal scons configuration to build openMalaria; it is not currently complete or platform-independant.
# TODO: boinc option
# TODO: print help
# TODO: separate build dir
# TODO: check for xsd or xsdcxx
# TODO: -O2 / -g, etc.

import os.path;
import subprocess;
import glob;
import dircache;

def which(program):
    import os
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def findSources(dPath):
  ret = []
  for l in dircache.listdir(dPath):
    d = dPath + '/' + l
    if os.path.isdir(d):
      ret = ret + findSources(d)
  ret = ret + glob.glob(dPath+'/*.cpp')
  ret = ret + glob.glob(dPath+'/*.C')
  ret = ret + glob.glob(dPath+'/*.cxx')
  return ret

def xsdFunc(target,source,env):
  outDir = os.path.dirname(str(target[0]))
  inFile = str(source[0])
  p = subprocess.Popen(env['XsdExec']+' cxx-tree --type-naming ucc --function-naming java --namespace-map =scnXml --generate-inline --generate-doxygen --output-dir '+outDir+' '+inFile, shell=True)
  return p.wait()

xsd = Builder(action = xsdFunc, src_suffix = '.xsd')
env = Environment(CPPPATH = ['include', 'xsdcxx'],
		  LIBPATH = ['lib'],
		  LIBS = ['gsl','gslcblas','xerces-c','z'])
env.Append(CPPSUFFIXES = [".d"])
env.Append(BUILDERS = {'XSD' : xsd})
#env.Append(CPPPATH = ['../boinc', '../boinc/api', '../boinc/lib'])
env.Append(CCFLAGS = ['-DWITHOUT_BOINC'])
conf = Configure(env)
# TODO: Checks for libraries, header files, etc. go here!
env['XsdExec'] = which('xsdcxx')
if env['XsdExec'] is None:
  env['XsdExec'] = which('xsd')
  if env['XsdExec'] is None:
    print 'Unable to find (code synthesis\') xsd program'
    Exit(1)

env = conf.Finish()

env.XSD(['xsdcxx/scenario.cxx','xsdcxx/scenario.hxx','xsdcxx/scenario.ixx'],'xsdcxx/scenario')
env.Program ('openMalaria', findSources('model') + findSources('xsdcxx'), PDB = 'openMalaria.pdb')
#env.Program ('openMalaria', ['model/gzstream.C'] + Glob('model/*.cpp') + Glob('model/*/*.cpp') + Glob('xsdcxx/*.c*'))

# Optimisations:

# check timestamp and if different checksum to decide whether to rebuild
Decider('MD5-timestamp')
# We can probably use this safely, if we need to.
#SetOption('implicit_cache', 1)
