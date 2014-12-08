#!/bin/env python2
# Note: ElementTree appears to have bugs in Python 3.3, but works in 2.7
# Most of this code works in both Python 2.x and 3.x, with a little compatibilty option:
from __future__ import print_function

# Tool for translating OpenMalaria XML files from one version to the next

import argparse, os, sys
import xml.etree.ElementTree as ET

LATEST_VERSION=33
KNOWN_VERSIONS=[32,33]    # all valid translation targets

class DocumentError(Exception):
    """Class for reporting errors in the input document."""
    def __init__(self,msg):
        Exception(self,msg)
        self.lines = []
    def __lshift__(self,line):
        self.lines.append(line)

class Translator:
    def __init__(self,f):
        self.tree = ET.parse(f)
        self.root = self.tree.getroot()
        self.schema_ver = int(self.root.get('schemaVersion', 0))
        if self.schema_ver not in KNOWN_VERSIONS:
            err = DocumentError('Unknown document or document version')
            err << 'Found version '+version
            err << 'Known versions: '+str(KNOWN_VERSIONS)
            raise err
        self.namespace = 'http://openmalaria.org/schema/scenario_' + str(self.schema_ver)
        if self.root.tag != self.tag_name():
            err=DocumentError('Unexpected tag name or namespace')
            err << 'Expected <'+self.tag_name()+'>'
            err << 'Found <'+self.root.tag+'>'
            raise err
    
    def write_to(self,f):
        ET.register_namespace('om',self.namespace)
        self.tree.write(f, encoding='UTF-8', xml_declaration=True)
    
    def update_to(self,target_version):
        while self.schema_ver < target_version:
            func_name = 'translate_' + str(self.schema_ver) + '_to_' + str(self.schema_ver+1)
            if not func_name in dir(self):      # do we have this function?
                raise NotImplementedError('Translation function '+func_name)
            getattr(self,func_name)()   # call it
            self.schema_ver += 1
            assert int(self.root.get('schemaVersion')) == self.schema_ver
    
    def tag_name(self):
        return '{' + self.namespace + '}scenario'
    
    def set_version_attrib(self,ver):
        oldVer = str(ver-1)
        newVer = str(ver)
        # This is the tag used by OpenMalaria to tell the version:
        self.root.set('schemaVersion', newVer)
        self.namespace = self.namespace.replace(oldVer, newVer)
        # This is the namespace 'om', which depends on the schema:
        self.root.tag = self.tag_name()
        # This is the location of that namespace. It may include the version
        # twice (in the namespace name and the schema file name).
        schemaLocName = '{http://www.w3.org/2001/XMLSchema-instance}schemaLocation'
        self.root.set(schemaLocName, self.root.get(schemaLocName).replace(oldVer, newVer))
    
    #def translate_32_to_33(self):
        #self.set_version_attrib(33)

def die(*args):
    print(*args, file=sys.stderr)
    sys.exit(1)

def readArgs():
    parser = argparse.ArgumentParser(description="This tool translates one or more "+
            "OpenMalaria scenario files (XML) from one version to the next.")
    parser.add_argument('files', metavar='FILE', nargs='*',
                        help="A file to translate. Multiple files may be specified.")   
    parser.add_argument('-t', '--target-version', type=int, action="store",
                        metavar='VER', default=LATEST_VERSION, choices=KNOWN_VERSIONS,
                        help="Target version for translation. Default: "+str(LATEST_VERSION))
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-d','--dest', metavar='DIR', action='store', default=None,
                        help='Destination directory for translated scenarios')
    group.add_argument('-i', '--in-place', action='store_true', default=False,
                        help='If set, files will be updated in-place (incompatible with -d option')
    args = parser.parse_args()
    if args.dest and not os.path.isdir(args.dest):
        die('Not a directory:', args.dest)
    return args

def main():
    args = readArgs()
    # args.target_version : target version number
    # args.files : list of files to translate
    # args.dest : destination directory or None
    # args.in_place : overwrite input files (True) or not (False)
    for f_path in args.files:
        try:
            if not os.path.isfile(f_path):
                die('Not a file:', f_path)
            with open(f_path, 'r') as f:
                translator = Translator(f)
            translator.update_to(args.target_version)
            if args.in_place:
                with open(f_path, 'w') as f:
                    translator.write_to(f)
            elif args.dest:
                path = os.path.join(args.dest, os.path.basename(f_path))
                with open(path, 'w') as f:
                    translator.write_to(f)
        except DocumentError,e:
            print('Error while translating',f_path)
            print(*e.args)
            for line in e.lines:
                print(line)

if __name__ == "__main__":
    main()
