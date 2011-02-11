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
import math
from readOutput import Keys, ValDict
from numbers import Number

measureNames = {
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

# for each combined measure, list included measures: id, name, colour
# encompasses all measures
combinedMeasures = [
    ('hosts',[(0,'all','black')]),
    ('infected hosts',[(1,'all','red'),(2,'expected','green'),(3,'patent','blue'),(43,'new','purple')]),
    ('sum log',[(4,'pyrogenic threshold','purple'),(5,'parasite density','green')]),
    ('total infections',[(6,'all','red'),(8,'patent','blue')]),
    ('transmitting humans',[(7,'sum(p(transmit))','green'),(26,'annual sum(p(transmit)) / annual EIR','blue')]),
    ('sum pyrog thres',[(10,'sum pyrogenic threshold','purple')]),
    ('treatments',[(11,'first line','red'),(12,'second line','purple'),(13,'hospital','orange')]),
    ('episodes',[(14,'UC','red'),(15,'severe','orange'),(27,'NMF','blue')]),
    ('sequelae',[(16,'all','red'),(24,'hospital','orange')]),
    ('deaths',[(17,'direct hospital','orange'),(18,'indirect','purple'),(19,'direct','red'),(41,'first day','brown'),(42,'first day hospital','green')]),
    ('intervention doses',[(20,'EPI vaccinations','red'),(22,'mass vaccinations','darkred'),(25,'IPT doses','brown')]),
    ('infant mortality rate',[(21,'IMR','red')]),
    ('hospital recoveries',[(23,'hospital recoveries','orange')]),
    ('inoculations',[(30,'all','red')]),
    ('feeding vectors',[(31,'emergence - N_v0','green'),(32,'all - N_v','blue'),(33,'infected - O_v','purple'),(34,'infectious - S_v','red')]),
    ('EIR (innocs/pers/year)',[(35,'requested','orange'),(36,'simulated','red')]),
    ('diagnostics used',[(39,'RDTs','purple'),(48,'microscopy','brown')]),
    ('drug usage (mg)',[(40,'oral','blue'),(49,'intrevenous','orange')]),
    ('vector interventions',[(44,'mass ITNs','darkgreen'),(45,'EPI ITNs','green'),(46,'mass IRS','darkred'),(47,'mass deterrents','brown')]),
    ('cohort delta',[(50,'added','orange'),(51,'removed','blue')])
]
def findMeasureGroup(m):
    for i in range(len(combinedMeasures)):
        for md in combinedMeasures[i][1]:
            if m==md[0]:
                return i
    raise KeyError("measure "+str(m)+" not in combinedMeasures")
def getMeasureLabel(mg,m):
    if mg is None:
        return measureNames[m]+" ("+str(m)+")"
    else:
        for md in combinedMeasures[mg][1]:
            if m==md[0]:
                return md[1]
        raise KeyError("measure "+str(m)+" not in combinedMeasures["+str(mg)+"]")
def getMeasureColour(mg,m):
    if mg is None:
        return 'blue'
    else:
        for md in combinedMeasures[mg][1]:
            if m==md[0]:
                return md[2]
        raise KeyError("measure "+str(m)+" not in combinedMeasures["+str(mg)+"]")

class MultiKey(object):
    __slots__ = ["mg","m","s","g","f"]
    def __init__(self,measureGroup=None,measure=None,survey=None,group=None,fileName=None):
        self.mg=measureGroup #int index in combinedMeasures
        self.m=measure #int
        self.s=survey #int
        self.g=group #string or int
        self.f=fileName #int
    def __and__(self,other):
        """Set each member from other is not set in self.
        Should complain if set in self and in other, but this is a bit unnecessary."""
        first=MultiKey(self.mg,self.m,self.s,self.g,self.f)
        if first.mg==None:
            first.mg=other.mg
        if first.m==None:
            first.m=other.m
        if first.s==None:
            first.s=other.s
        if first.g==None:
            first.g=other.g
        if first.f==None:
            first.f=other.f
        return first
    @staticmethod
    def fromMeasureGroups(measureGroups):
        r= [MultiKey(measureGroup=measureGroup) for measureGroup in measureGroups]
        return r
    @staticmethod
    def fromMeasures(measures):
        r= [MultiKey(measure=measure) for measure in measures]
        return r
    @staticmethod
    def fromSurveys(surveys):
        return [MultiKey(survey=survey) for survey in surveys]
    @staticmethod
    def fromGroups(groups):
        return [MultiKey(group=group) for group in groups]
    @staticmethod
    def fromFiles(files):
        return [MultiKey(fileName=fl) for fl in files]
    @staticmethod
    def expand(r,x):
        lx=len(x)
        assert lx>0, "expected len>0"
        r *= lx # shallow copy: all references point to same object
        for i in range(len(r)):
            r[i] = r[i] & x[i%lx]
    def __str__(self):
        return self.label(MultiKey(),True,None)
    def label(self,base,replaceFN,valDict):
        r=""
        if self.m!=base.m:
            r=getMeasureLabel(self.mg,self.m)
        if self.s!=base.s:
            if len(r):
                r+=","
            r+="survey "+str(self.s)
        if self.g!=base.g:
            if len(r):
                r+=","
            r+="group "+str(self.g)
        if self.f!=base.f:
            if len(r):
                r+=","
            if replaceFN:
                r+="run "+str(self.f)
            else:
                r+=valDict.getFileName(self.f)
        return r

class Plotter(object):
    def __init__(self,keys):
        self.values=ValDict(keys)
    def read(self,fileName,measures):
        self.values.read(fileName,measures)
    def plot(self,am,replaceFN,showLegends,s,g,f):
        x_axis=Keys.NONE
        x_label=""
        if s=="x-axis":
            assert x_axis==Keys.NONE, "dual assignment to x-axis!"
            x_axis=Keys.SURVEY
            x_label="survey"
        if g=="x-axis":
            assert x_axis==Keys.NONE, "dual assignment to x-axis!"
            x_axis=Keys.GROUP
            x_label="group"
        if f=="x-axis":
            assert x_axis==Keys.NONE, "dual assignment to x-axis!"
            x_axis=Keys.FILE
            x_label="file"
            x=self.values.getFileNames(replaceFN)
        assert x_axis!=Keys.NONE, "nothing to plot on x-axis!"
        
        lines=set()
        if s=="line": lines.add(Keys.SURVEY)
        if g=="line": lines.add(Keys.GROUP)
        if f=="line": lines.add(Keys.FILE)
        
        plots=[MultiKey()]
        if am:
            measureGroups=dict()
            for m in self.values.getMeasures():
                mg=findMeasureGroup(m)
                if mg in measureGroups:
                    measureGroups[mg].append(m)
                else:
                    measureGroups[mg]=list([m])
            #note: here "measure" is used to represent "measure group"
            MultiKey.expand(plots, MultiKey.fromMeasureGroups(measureGroups.keys()))
        else:
            MultiKey.expand(plots, MultiKey.fromMeasures(self.values.getMeasures()))
        if s=="plot":
            MultiKey.expand(plots, MultiKey.fromSurveys(self.values.getSurveys()))
        if g=="plot":
            raise Exception("Cannot separate group by plot!")
        if f=="plot":
            MultiKey.expand(plots, MultiKey.fromFiles(self.values.getFiles()))
        
        n=len(plots)
        d1=int(math.ceil(math.sqrt(float(n))))
        d2=int(math.ceil(float(n)/float(d1)))
        
        fig = plt.figure(1)
        i=1
        for plot in plots:
            subplot = fig.add_subplot(d1,d2,i)
            i+=1
            
            m=plot.m
            if am:
                m=measureGroups[plot.mg][0] # take first — not necessary correct but we need a measure
            
            if x_axis==Keys.GROUP:
                x=self.values.getGroups(m)
                x_label=self.values.getGroupLabel(m)
            elif x_axis==Keys.SURVEY:
                x=self.values.getSurveys(m)
            #else x was set previously
            
            subplot.set_title(plot.label(MultiKey(measure=plot.m),replaceFN,self.values))
            # Don't set x-label: is always the same and gets in the way of title
            #subplot.set_xlabel(x_label)
            if am:
                subplot.set_ylabel(combinedMeasures[plot.mg][0])
                
                pLines=[plot]
                measures=measureGroups[plot.mg]
                MultiKey.expand(pLines, MultiKey.fromMeasures(measures))
            else:
                subplot.set_ylabel(getMeasureLabel(plot.mg,plot.m))
                
                pLines=[plot]
            
            #TODO: colour is not unique for any key other than measures from measure groups
            # i.e. we need to somehow manipulate colours to make them unique
            if Keys.SURVEY in lines:
                MultiKey.expand(pLines, MultiKey.fromSurveys(self.values.getSurveys(m)))
            if Keys.GROUP in lines:
                MultiKey.expand(pLines, MultiKey.fromGroups(self.values.getGroups(m)))
            if Keys.FILE in lines:
                MultiKey.expand(pLines, MultiKey.fromFiles(self.values.getFiles()))
            
            if len(x)>1 and isinstance(x[0], Number): # draw an xy line chart
                plotted=list()
                for pLine in pLines:
                    xKeys=[pLine]
                    if x_axis==Keys.SURVEY:
                        MultiKey.expand(xKeys, MultiKey.fromSurveys(self.values.getSurveys(m)))
                    elif x_axis==Keys.GROUP:
                        MultiKey.expand(xKeys, MultiKey.fromGroups(self.values.getGroups(m)))
                    elif x_axis==Keys.FILE:
                        MultiKey.expand(xKeys, MultiKey.fromFiles(self.values.getFiles()))
                    y=[self.values.get(k.m,k.s,k.g,k.f) for k in xKeys]
                    colour=getMeasureColour(pLine.mg,pLine.m)
                    try:
                        plotted.append(subplot.plot(x,y,colour))
                    except ValueError,e:
                        print "Bad plot values (script error):"
                        print "x:",x
                        print "y:",y
                
                if showLegends and len(plotted)>1:
                    legends=[pLine.label(plot,replaceFN,self.values) for pLine in pLines]
                    subplot.legend(plotted,legends,'upper right')
            else: #one x-coord or non-numeric x-coords: draw a bar chart
                plotted=list()
                xind=numpy.arange(len(x))
                width=0.9/len(pLines)
                xincr=0.0
                for pLine in pLines:
                    xKeys=[pLine]
                    if x_axis==Keys.SURVEY:
                        MultiKey.expand(xKeys, MultiKey.fromSurveys(self.values.getSurveys(m)))
                    elif x_axis==Keys.GROUP:
                        MultiKey.expand(xKeys, MultiKey.fromGroups(self.values.getGroups(m)))
                    elif x_axis==Keys.FILE:
                        MultiKey.expand(xKeys, MultiKey.fromFiles(self.values.getFiles()))
                    y=[self.values.get(k.m,k.s,k.g,k.f) for k in xKeys]
                    colour=getMeasureColour(pLine.mg,pLine.m)
                    try:
                        plotted.append(subplot.bar(xind+xincr,y,width,color=colour))
                    except ValueError,e:
                        print "Bad plot values (script error):"
                        print "x:",xind+xincr
                        print "y:",y
                    xincr+=width
                
                if showLegends and len(plotted)>1:
                    plots=[p[0] for p in plotted]
                    legends=[pLine.label(plot,replaceFN,self.values) for pLine in pLines]
                    subplot.legend(plots,legends,'upper right')
                if len(x)>1:
                    subplot.set_xticks(xind+0.45)
                    subplot.set_xticklabels(x)
                subplot.set_xlim((-0.1,len(x)))
        
        plt.show()

def main(args):
    parser = OptionParser(usage="Usage: %prog [options] FILES",
            description="""Plots results from an OpenMalaria (surveys) output
file by time. Currently no support for simultaeneously handling
multiple files or plotting according to age group.

Valid targets for plotting keys are: none (key is aggregated), x-axis, plot, line.
If no key is set to the x-axis, the first unassigned of survey, group, file will be
assigned to the x-axis.""",version="%prog 0.1")
    
    parser.add_option("-m","--measures", action="store", type="string", dest="measures", default="",
            help="Plot only MEASURES (comma-separated list of numbers)")
    parser.add_option("-a","--no-auto-measures", action="store_false", dest="am", default=True,
            help="Don't automatically put similar measures on the same plot.")
    #parser.add_option("-m","--measure", action="store", type="choice", dest="m", default="plot",
            #choices=["none","x-axis","plot","line"],help="How to plot measures")
    parser.add_option("-s","--survey", action="store", type="choice", dest="s", default="none",
            choices=["none","x-axis","plot","line"],help="How to plot surveys")
    parser.add_option("-g","--group", action="store", type="choice", dest="g", default="none",
            choices=["none","x-axis","plot","line"],help="How to plot (age) groups")
    parser.add_option("-f","--file", action="store", type="choice", dest="f", default="plot",
            choices=["none","x-axis","plot","line"],help="How to plot outputs from different files")
    parser.add_option("-n","--file-names", action="store_true", dest="fn", default=False,
            help="Use full file names instead of replacing with short versions")
    parser.add_option("-l","--no-legends", action="store_true", dest="nl", default=False,
            help="Turn off legends (sometimes they hide parts of the plots)")
    
    (options, others) = parser.parse_args(args=args[1:])
    if len(others)==0:
        parser.print_usage()
        return 1
    
    measures=options.measures.split(",")
    if(measures[-1].strip()==""):
        measures=measures[:-1]
    measures=[int(m) for m in measures]
    
    if options.s != "x-axis" and options.g != "x-axis" and options.f != "x-axis":
        if options.s == "none":
            options.s="x-axis"
        elif options.g == "none":
            options.g="x-axis"
        elif options.f == "none":
            options.f="x-axis"
        else:
            print "Error: nothing assigned to x-axis!"
            return 1
    
    keys=set()
    keys.add(Keys.MEASURE)
    if options.s != "none": keys.add(Keys.SURVEY)
    if options.g != "none": keys.add(Keys.GROUP)
    if options.f != "none": keys.add(Keys.FILE)
    
    plotter=Plotter(keys)
    
    for output in others:
        plotter.read(output,measures)
    
    plotter.plot(options.am,not options.fn,not options.nl,options.s,options.g,options.f)
    
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
