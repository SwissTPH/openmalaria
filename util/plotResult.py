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

from optparse import OptionParser
import numpy
import matplotlib.pyplot as plt
from matplotlib.colors import cnames
import sys
import string
import math as m

measure_names = {
    0 : 'nHost',
    1 : 'nInfect',
    2 : 'nExpectd',
    3 : 'nPatent',
    4 : 'sumLogPyrogenThres',
    5 : 'sumlogDens',
    6 : 'totalInfs',
    7 : 'nTransmit',
    8 : 'totalPatentInf',
    9 : 'contrib',
    10 : 'sumPyrogenThresh',
    11 : 'nTreatments1',
    12 : 'nTreatments2',
    13 : 'nTreatments3',
    14 : 'nUncomp',
    15 : 'nSevere',
    16 : 'nSeq',
    17 : 'nHospitalDeaths',
    18 : 'nIndDeaths',
    19 : 'nDirDeaths',
    20 : 'nEPIVaccinations',
    21 : 'allCauseIMR',
    22 : 'nMassVaccinations',
    23 : 'nHospitalRecovs',
    24 : 'nHospitalSeqs',
    25 : 'nIPTDoses',
    26 : 'annAvgK',
    27 : 'nNMFever',
    30 : 'innoculationsPerAgeGroup',
    28 : 'innoculationsPerDayOfYear',
    29 : 'kappaPerDayOfYear',
    31 : 'Vector_Nv0',
    32 : 'Vector_Nv',
    33 : 'Vector_Ov',
    34 : 'Vector_Sv',
    35 : 'Vector_EIR_Input',
    36 : 'Vector_EIR_Simulated',
    39 : 'Clinical_RDTs',
    40 : 'Clinical_DrugUsage',
    41 : 'Clinical_FirstDayDeaths',
    42 : 'Clinical_HospitalFirstDayDeaths',
    43 : 'nNewInfections',
    44 : 'nMassITNs',
    45 : 'nEPI_ITNs',
    46 : 'nMassIRS',
    47 : 'nMassVA',
    48 : 'Clinical_Microscopy',
    49 : 'Clinical_DrugUsageIV',
    50 : 'nAddedToCohort',
    51 : 'nRemovedFromCohort',
    52 : 'nMDAs',
    53 : 'nNmfDeaths',
    54 : 'nAntibioticTreatments',
    }

class Value(object):
    def __init__(self):
        self.n=0
        self.v=0.0
    
    def add(self,val):
        self.v+=val
        self.n+=1
    
    def get(self):
        if self.n > 0:
            return self.v
            #return self.v/self.n
        else:
            return 1e10000 * 0 # NaN

class SimOutputs(object):
    def __init__(self,outFile):
        measures=set()
        maxSurvey=0
        f=open(outFile)
        for line in f:
            items=string.split(line,"\t")
            if (len(items) != 4):
                raise Exception("expected 4 items on line: "+line)
            measures.add(int(items[2]))
            maxSurvey=max(maxSurvey,int(items[0]))
        self.m=dict()
        for measure in measures:
            self.m[measure]=[Value() for i in range(maxSurvey+1)]
        f.close()
        f=open(outFile)
        for line in f:
            items=string.split(line,"\t")
            if (len(items) != 4):
                raise Exception("expected 4 items on line")
            self.m[int(items[2])][int(items[0])].add(float(items[3]))
    
    def plot(self,measure,subplot):
        x=range(len(self.m[measure]))
        y=[v.get() for v in self.m[measure]]
        #different colours (possibly also white): color=cnames.values()[measure]
        p=subplot.plot(x,y)

def plot(outFile,measures):
    output=SimOutputs(outFile)
    if(len(measures)==0):
        measures=output.m.keys()
        measures.sort()
    
    fig = plt.figure(1)
    n=len(measures)
    d1=int(m.ceil(m.sqrt(float(n))))
    d2=int(m.ceil(float(n)/float(d1)))
    i=1
    for measure in measures:
        p = fig.add_subplot(d1,d2,i)
        p.set_xlabel('survey')
        try:
            p.set_ylabel(measure_names[measure]+' ('+str(measure)+')')
        except KeyError:
            p.set_ylabel('measure '+str(measure))
        output.plot(measure,p)
        i+=1
    plt.show()

def main(args):
    parser = OptionParser(usage="Usage: %prog [options] FILE",
            description="""Plots results from an OpenMalaria (surveys) output
file by time. Currently no support for simultaeneously handling
multiple files or plotting according to age group.""",version="%prog 0.1")
    
    parser.add_option("-m", action="store", type="string", dest="measures",
            help="Plot only measures X (comma-separated list of numbers")
    
    (options, others) = parser.parse_args(args=args)
    if len(others)==0:
        parser.print_usage()
        return 1
    
    options.ensure_value("measures", "")
    measures=options.measures.split(",")
    if(measures[-1].strip()==""):
        measures=measures[:-1]
    
    for output in others[1:]:
        plot(output,[int(m) for m in measures])
    
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
