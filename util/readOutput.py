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

import string
import unittest

class Keys:
    NONE=0
    MEASURE=1
    SURVEY=2
    GROUP=3
    FILE=4
    
    all=set([NONE,MEASURE,SURVEY,GROUP,FILE])
    
    def fromString(str):
        if str=="none":
            return NONE
        elif str=="measure":
            return MEASURE
        elif str=="survey":
            return SURVEY
        elif str=="group":
            return GROUP
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
        self.assert_ (self.a1.a == self.a2.a)
        self.assert_ (self.a1 == self.a2)
        self.assert_ (self.b1 == self.b2)
        self.assert_ (self.a1 != self.a3)
    def testHash (self):
        self.assert_ (self.a1.__hash__() == self.a2.__hash__())
        self.assert_ (self.a1.__hash__() != self.a3.__hash__()) # actually, hash collisions are possible
        self.assert_ (self.b1.__hash__() == self.b2.__hash__())

def isAgeGroup(measure):
    if measure in set([7,21,26,31,32,33,34,35,36,39,40,48,49]):
        return False
    return True

class MeasureDict(object):
    def add(self,survey,group,f,value):
        """survey:int, group:string or int, f:int, value:float"""
        pass
    def get(self,survey,group,f):
        """survey:int, group:string, f:string"""
        pass
    def getGroups(self):
        pass
class MeasureAGDict(object):
    def __init__(self):
        self.nAGroups = 0
        #(list by fileID) of (list by survey) of (list by group)
        self.v=list()
        self.groupLabel="age group"
    def add(self,survey,group,f,value):
        """survey: int, group: string, f: int"""
        ag = int(group) # force type as int
        if ag >= self.nAGroups:
            self.nAGroups = ag+1
        
        while f>=len(self.v):
            self.v.append(list())
        surveys=self.v[f]
        while survey >= len(surveys):
            surveys.append(list())
        groups=surveys[survey]
        while ag >= len(groups):
            groups.append(0.0)
        groups[ag] += value
    def get(self,survey,group,f):
        try:
            return self.v[f][survey][int(group)]
        except KeyError:
            return 1e1000 - 1e0000 # NaN
    def getGroups(self):
        return range(0,self.nAGroups)
class MeasureOGDict(object):
    def __init__(self,m):
        self.groups=set() # record of all groups this measure has
        #(list by fileID) of (list by survey) of (dict by group)
        self.v=list()
        if m in set([40,49]):
            self.groupLabel="drug ID"
        elif m>=31 and m<=34:
            self.groupLabel="vector species"
        else:
            self.groupLabel="(none)"
    def add(self,survey,group,f,value):
        """survey: int, group: string, f: int"""
        group=str(group) # force same type to make sure hash sums match
        self.groups.add(group)
        
        while f>=len(self.v):
            self.v.append(list())
        surveys=self.v[f]
        while survey >= len(surveys):
            surveys.append(dict())
        groups=surveys[survey]
        groups[group] = groups.get(group,0.0) + value
    def get(self,survey,group,f):
        group=str(group)
        try:
            return self.v[f][survey][group]
        except KeyError:
            print "can't find:",f,survey,group
            print "have:",self.v
            return 1e1000 - 1e1000 # NaN
    def getGroups(self):
        return list(self.groups)

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
    
    def read(self,fileName,measures):
        """Read from fileName. If measures is non-empty, only read these measures."""
        if len(measures)==0:
            measures=None
        kS = Keys.SURVEY not in self.aggregateKeys
        kG = Keys.GROUP not in self.aggregateKeys
        if Keys.FILE not in self.aggregateKeys:
            assert fileName not in self.files, "Reading same file twice?"
            fID = len(self.files)
            self.files.append(fileName)
        else:
            fID = 0
        s = 0
        g = 0
        fileObj = open(fileName, 'r')
        nErrs=0
        for line in fileObj:
            items=string.split(line)
            if (len(items) != 4):
                print "expected 4 items on line; found (following line):"
                print line
                nErrs+=1
                if nErrs>5:
                    raise Exception ("Too many errors reading "+fileName)
                continue
            
            m=int(items[2])
            if measures!=None and m not in measures:
                continue
            self.measures.add(m)
            i=len(self.values)
            while m >= i:
                self.values.append(MeasureAGDict() if isAgeGroup(i) else MeasureOGDict(m))
                i+=1
            if kS:
                s=int(items[0])
                self.nSurveys=max(self.nSurveys,s)
            if kG:
                g=items[1]
            self.values[m].add(s,g,fID,robustFloat(items[3]))
    
    def getFiles(self):
        return range(len(self.files))
    def getFileName(self,n):
        return self.files[n]
    def getFileNames(self,replaceFN):
        if replaceFN:
            return ["run "+str(n) for n in range(len(self.files))]
        else:
            return self.files
    def getMeasures(self):
        return list(self.measures)
    def getSurveys(self,m):
        """takes measure no. This does something special for measure 21, otherwise
        normal behaviour, so passing a const like 0 is fine."""
        if m==21:
            return [1]
        else:
            return range(1,self.nSurveys+1)
    def getGroups(self,measure):
        return self.values[measure].getGroups()
    def getGroupLabel(self,measure):
        return self.values[measure].groupLabel
    def get(self,m,s,g,f):
        if s==None:
            s=0
        if g==None:
            g=0
        if f==None:
            f=0
        return self.values[m].get(s,g,f)


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
    """Return a dict of entries read from file. Keys have type Multi3Keys.
    
    Note: ValDict is probably more efficient due to use of arrays over dicts."""
    values=dict()
    fileObj = open(fname, 'r')
    for line in fileObj:
        items=string.split(line)
        if (len(items) != 4):
            print "expected 4 items on line; found (following line):"
            print line
            continue
            
        key=Multi3Keys(int(items[0]),items[1],int(items[2]))
        values[key]=robustFloat(items[3])
    return values

if __name__ == '__main__':
    unittest.main()
