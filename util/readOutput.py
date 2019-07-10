#!/usr/bin/env python3
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

import unittest

class Keys:
    NONE=0
    MEASURE=1
    SURVEY=2
    GROUP=3
    COHORT=4
    GENOTYPE=5
    FILE=6
    
    all=set([NONE,MEASURE,SURVEY,GROUP,COHORT,GENOTYPE,FILE])
    
    def fromString(str):
        if str=="none":
            return NONE
        elif str=="measure":
            return MEASURE
        elif str=="survey":
            return SURVEY
        elif str=="group":
            return GROUP
        elif str=="cohort":
            return COHORT
        elif str=="genotype":
            return GENOTYPE
        elif str=="file":
            return FILE
        else:
            raise Exception("invalid key: "+str)

class Multi3Keys(object):
    __slots__=["a","b","c"]
    """Class combining three keys into a single key."""
    def __init__(self,a,b,c):
        self.a = a
        self.b = b
        self.c = c
    def __eq__(self,other):
        return (self.a == other.a) and (self.b == other.b) and (self.c == other.c)
    def __hash__(self):
        return self.a.__hash__() ^ self.b.__hash__() ^ self.c.__hash__()

class TestMultiKeys (unittest.TestCase):
    def setUp(self):
        self.a1 = Multi3Keys(2,4,0);
        self.a2 = Multi3Keys(2,4,0);
        self.a3 = Multi3Keys(2,2,0);
        self.b1 = Multi3Keys(0,"abc",5);
        self.b2 = Multi3Keys(0,"abc",5);
    def testEq (self):
        self.assertTrue (self.a1.a == self.a2.a)
        self.assertTrue (self.a1 == self.a2)
        self.assertTrue (self.b1 == self.b2)
        self.assertTrue (self.a1 != self.a3)
    def testHash (self):
        self.assertTrue (self.a1.__hash__() == self.a2.__hash__())
        self.assertTrue (self.a1.__hash__() != self.a3.__hash__()) # actually, hash collisions are possible
        self.assertTrue (self.b1.__hash__() == self.b2.__hash__())

def isAgeGroup(measure):
    if measure in set([7,9,21,25,26,28,29,31,32,33,34,35,36,39,40,47,48,49,50,51,54]):
        return False
    return True

class MeasureDict(object):
    def __init__(self,m):
        self.nGroups = 0
        self.nCohorts = 0
        self.nGenotypes = 0
        #(list by fileID) of (list by survey) of (list by group) of (list by cohort) of (list by genotype)
        self.v=list()
        if m>=31 and m<=34:
            self.groupLabel="vector species"
        else:
            self.groupLabel="age group"
    def add(self,survey,group,cohort,genotype,f,value):
        """survey:int, group:int, cohort:int, genotype:int, f:int, value:float"""
        self.nGroups = max(self.nGroups,group+1)
        self.nCohorts = max(self.nCohorts,cohort+1)
        self.nGenotypes = max(self.nGenotypes,genotype+1)
        
        while f>=len(self.v):
            self.v.append(list())
        surveys=self.v[f]
        while survey >= len(surveys):
            surveys.append(list())
        groups=surveys[survey]
        while group >= len(groups):
            groups.append(list())
        cohorts=groups[group]
        while cohort >= len(cohorts):
            cohorts.append(list())
        genotypes=cohorts[cohort]
        while genotype >= len(genotypes):
            genotypes.append(0.0)
        genotypes[genotype] += value
    def get(self,survey,group,cohort,genotype,f):
        try:
            return self.v[f][survey][group][cohort][genotype]
        except LookupError:
            return 1e1000 - 1e0000 # NaN
    def getGroups(self):
        return list(range(0,self.nGroups))
    def getCohorts(self):
        return list(range(0,self.nCohorts))
    def getGenotypes(self):
        return list(range(0,self.nGenotypes))

def stringIndexAllMatch(strs,ind,char):
    for s in strs:
        if s[ind] != char:
            return False
    return True

class ValDict (object):
    """Class looking like a dictionary of outputs, but supporting aggregation
    and keeping lists of all keys.
    
    types of keys (see Keys) to separate by (others are aggregated): 
    """
    def __init__(self,keys):
        self.aggregateKeys = Keys.all - keys
        self.nSurveys=0 # set to max survey number; indecies are +1
        self.values=list() #key: measure number
        self.measures=set() #set of used measures
        self.files=list()
    
    def read(self,fileName,filterExpr,exprDebug):
        """Read from fileName. If measures is non-empty, only read these measures."""
        def filterFun(f,m,s,g,c,gt):
            r=eval(filterExpr)
            if exprDebug:
                print(("f="+str(f),"m="+str(m),"s="+str(s),"g="+str(g),"c="+str(c),"g="+str(g)+":",r))
            return r
        aKS = Keys.SURVEY in self.aggregateKeys
        aKG = Keys.GROUP in self.aggregateKeys
        aKC = Keys.COHORT in self.aggregateKeys
        aKGT = Keys.GENOTYPE in self.aggregateKeys
        if Keys.FILE not in self.aggregateKeys:
            assert fileName not in self.files, "Reading same file twice?"
            fID = len(self.files)
            self.files.append(fileName)
        else:
            fID = 0
        fileObj = open(fileName, 'r')
        nErrs=0
        for line in fileObj:
            items=line.split()
            if (len(items) != 4):
                print("expected 4 items on line; found (following line):")
                print(line)
                nErrs+=1
                if nErrs>5:
                    raise Exception ("Too many errors reading "+fileName)
                continue
            
            m=int(items[2])
            s=int(items[0])
            g=int(items[1])
            gt = g / 1000000 # genotype
            g = g - 1000000*gt
            c = g / 1000   # cohort
            g = g - 1000*c
            if not filterFun(fileName,m,s,g,c,gt):
                continue
            if aKS:
                s=0
            if aKG:
                g=0
            if aKC:
                c=0
            if aKGT:
                gt=0
            i=len(self.values)
            while m >= i:
                self.values.append(MeasureDict(m))
                i+=1
            self.measures.add(m)
            self.nSurveys=max(self.nSurveys,s)
            self.values[m].add(s,g,c,gt,fID,robustFloat(items[3]))
    
    def getFiles(self):
        return list(range(len(self.files)))
    def getFileName(self,n):
        return self.files[n]
    def getFileNames(self,replaceFN):
        if not hasattr(self,"fnIndex"):
            if len(self.files) <= 1:
                self.fnIndex=0
            else:
                i=0
                try:
                    while True:
                        c=self.files[0][i]
                        if not stringIndexAllMatch(self.files[1:],i,c):
                            break
                        i+=1
                except IndexError:
                    pass
                self.fnIndex=i
        longNames=[f[self.fnIndex:] for f in self.files]
        if replaceFN and max([len(n) for n in longNames])>8:
            return ["file "+str(n+1) for n in range(len(self.files))]
        else:
            return longNames
    def getMeasures(self):
        return list(self.measures)
    def getSurveys(self,m):
        """takes measure no. This does something special for measure 21, otherwise
        normal behaviour, so passing a const like 0 is fine."""
        if m==21:
            return [1]
        else:
            return list(range(1,self.nSurveys+1))
    def getAllGroups(self):
        groups=set()
        for x in self.values:
            for y in x.getGroups():
                groups.add(y)
        return groups
    def getGroups(self,measure):
        return self.values[measure].getGroups()
    def getCohorts(self,measure):
        return self.values[measure].getCohorts()
    def getGenotypes(self,measure):
        return self.values[measure].getGenotypes()
    def getGroupLabel(self,measure):
        return self.values[measure].groupLabel
    def get(self,m,s,g,c,gt,f):
        if s==None:
            s=0
        if g==None:
            g=0
        if c==None:
            c=0
        if gt==None:
            gt=0
        if f==None:
            f=0
        return self.values[m].get(s,g,c,gt,f)


#http://stackoverflow.com/questions/2974124/reading-floating-point-numbers-with-1-qnan-values-in-python
def robustFloat(s):
    """Return an NaN instead of throwing."""
    try:
        return float(s)
    except ValueError:
        if 'nan' in s.lower():
            return 1e1000-1e1000 # NaN
        else:
            raise

def readEntries (fname):
    """Return a dict of entries read from file. Keys have type Multi3Keys,
    where a corresponds to measure, b to survey and c to group.
    
    Note: ValDict is probably more efficient due to use of arrays over dicts."""
    values=dict()
    fileObj = open(fname, 'r')
    for line in fileObj:
        items=line.split()
        if (len(items) != 4):
            print("expected 4 items on line; found (following line):")
            print(line)
            continue
            
        key=Multi3Keys(int(items[2]),int(items[0]),int(items[1]))
        values[key]=robustFloat(items[3])
    return values

if __name__ == '__main__':
    unittest.main()
