#!/usr/bin/python
# -*- coding: utf-8 -*-

# This file is part of OpenMalaria.
# 
# Copyright (C) 2005-2011 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

# This is intended to be a generalization of plotResults.py which can use
# combined outputs as produced by resultAssimilator.py. It is currently
# very incomplete.

import string
import os.path
import sys
from optparse import OptionParser

class DataProcessor(object):
    def __init__(self,sep):
        self.sep = sep
        self.plot=dict()
    def add(p,l,x,v):
        self.plot[p][l][x]
    def read(self,file_name):
        f=open(file_name)
        line=f.readline().strip().split(self.sep)
        assert line[-1] == 'value'
        

def main(args):
    parser = OptionParser(usage="Usage: %prog [options] FILE",
            description="""FILE is assumed to be a csv-style file. The first line should be a header, the last entry of which should be 'value'.""",version="%prog 0.1")
    
    parser.add_option("-e","--filter", action="store", type="string", dest="filterExpr", default="m!=0",
            help="Filter entries read according to this expression (i.e. values are included when this returns true). Parameters available are the column names (from header).. Examples: 'True', 'm!=0' (default), 'm in [11,12,13]', 's > 73 and m!=0'.")
    parser.add_option("--debug-filter", action="store_true", dest="debugFilter", default=False,
            help="Each time FILTEREXPR is called, print input values and output. Warning: will print a lot of data!")
    parser.add_option("-s","--separator", action="store", type="string", dest="sep", default="\t",
            help="Specify the column separator (defaults to \\t).")
    
    (options, others) = parser.parse_args(args=args[1:])
    if len(others)!=1:
        parser.print_usage()
        return 1
    
    dp = DataProcessor(options.sep,others[0])
    dp.
    
    return 0
def main(args):
    parser = OptionParser(usage="Usage: %prog [options] SCENARIOS.CSV RESULTS_DIR OUTPUT_FILE",
            description="""""",version="%prog 0.1")
    
    (options, others) = parser.parse_args(args=args[1:])
    if len(others)!=3 or not os.path.isdir(others[1]):
        parser.print_usage()
        return 1
    
    sweeps,swArmIds = readSwArmIds(others[0])
    
    outFile=open(others[2],'w')
    resultDir=others[1]
    
    for sweep in sweeps:
        outFile.write(sweep)
        outFile.write('\t')
    outFile.write('sumNewInfs\n')
    
    for f,k in swArmIds.iteritems():
        lineStart=''
        for x in k:
            lineStart+=x+'\t'
        resPath=os.path.join(resultDir,f.replace('.xml','.txt'))
        if os.path.isfile(resPath):
            res=open(resPath)
            newInfs=0.
            for line in res:
                lineS=line.strip().split('\t')
                if int(lineS[0])>2 and int(lineS[2])==43:
                    newInfs+=float(lineS[3])
            outFile.write(lineStart)
            outFile.write(str(newInfs)+'\n')
            res.close()
        else:
            pass # missing result
    
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
