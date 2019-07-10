#!/usr/bin/env python3
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

import os.path
import sys
from optparse import OptionParser

def readSwArmIds(fileName):
    """Read scenarios.csv file of all arm ids for each sweep per scenario."""
    handle=open(fileName)
    
    # get headers, skipping file name
    headers=handle.readline().strip().split(',')[1:]
    ret=dict()
    for line in handle:
        csi=line.strip().split(',')
        ret[csi[0]] = csi[1:]
    return headers,ret

def main(args):
    parser = OptionParser(usage="Usage: %prog [options] SCENARIOS.CSV RESULTS_DIR OUTPUT_FILE",
            description="""Given a .csv file associating scenario file names
            with DB scenario IDs (sce_id), SCENARIOS.CSV, and a directory of results,
            RESULTS_DIR, where results have the name of the scenario file but
            with .xml substituted for .txt, a tab-separated file of combined
            outputs, OUTPUT_FILE, is written, where each line contains sce_id,
            survey, group and measure identifiers, and a value.""",version="%prog 0.1")
    
    (options, others) = parser.parse_args(args=args[1:])
    if len(others)!=3 or not os.path.isdir(others[1]):
        parser.print_usage()
        return 1
    
    sweeps,swArmIds = readSwArmIds(others[0])
    sce_id_ind = next((i for i in range(len(sweeps)) if sweeps[i] == 'sce_id'), None)
    if sce_id_ind is None:
        raise Exception("unable to find column for sce_id!")
    # sce_id_ind is zero-indexed and misses first column, so add 2:
    print(("Column for sce_id:",sce_id_ind+2))
    
    outFile=open(others[2],'w')
    resultDir=others[1]
    
    for f,k in list(swArmIds.items()):
        sce_id=k[sce_id_ind]
        lineStart=sce_id+'\t'
        resPath=os.path.join(resultDir,f+'.txt')
        if os.path.isfile(resPath):
            res=open(resPath)
            for line in res:
                outFile.write(lineStart)
                outFile.write(line)
            res.close()
        else:
            pass # missing result
    
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
