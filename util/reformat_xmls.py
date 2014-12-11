#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This file is part of OpenMalaria.

Copyright (C) 2005-2014 Swiss Tropical Institute

OpenMalaria is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
"""
"""
This script is used to reformat the xml files if you're experiencing problems
with whitespaces. You will need to have PyXML installed.
Usage: ./reformat_xmls.py -p folder_path | file_path
"""

from xml.dom.minidom import parse, parseString
# xml.dom.ext comes from PyXML which is no longer maintained an not always available :-(
from xml.dom.ext import PrettyPrint
from os.path import exists
from os.path import isdir
from os.path import isfile
import os
import sys
import getopt

def write_to_file(doc, output_file_path):
    file_object = open(output_file_path, "w")
    PrettyPrint(doc, file_object)
    file_object.close()

def parse_and_write(input_file_path):
    basename, extension = input_file_path.split('.')
    if(extension == 'xml'):
        doc = parse(input_file_path)
        write_to_file(doc, input_file_path)


opts, extraparams = getopt.getopt(sys.argv[1:], 'p', 'path')

for o, p in opts:
    if o in ['-p', '--path']:
        path = extraparams[0]

if exists(path) :
    if isdir(path):
        dirList = os.listdir(path)
        for fname in dirList:
            parse_and_write(path + "/" + fname)
    elif isfile(path):
            parse_and_write(path)
else:
    print("Please enter a valid file or folder path...")
