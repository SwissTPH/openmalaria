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

import sys
import string
from optparse import OptionParser
from approxEqual import ApproxSame

REL_PRECISION=1e-6
ABS_PRECISION=1e-6


def charEqual (fn1,fn2):
    MAX=10*1024
    f1 = open(fn1,'r')
    f2 = open(fn2,'r')
    while True:
        s1 = f1.read(MAX)
        s2 = f2.read(MAX)
        if (len(s1)==0) or (len(s2)==0):
            # end of one or both files; equal if it's the end of both
            return len(s1) == len(s2)
        if s1 != s2:
            return False

def readSums(fn):
    """return list of titles and list of column sums"""
    f=open(fn,'r')
    titles=[s.strip() for s in f.readline().split('\t')]
    if titles==['##','##']:
        # first line is auxilliary header...
        titles=[s.strip() for s in f.readline().split('\t')]
    cols=[0.0 for x in range(len(titles))]
    for line in f:
        i=0
        for item in line.split('\t'):
            cols[i] += float(item.strip())
            i+=1
    return titles,cols

def sumsApproxEq (fn1,fn2):
    """ return true if approx equal, false otherwise
    also prints about unequal stuff"""
    t1,vec1=readSums(fn1)
    t2,vec2=readSums(fn2)
    approxSame = ApproxSame(REL_PRECISION,ABS_PRECISION)
    if t1 != t2:
        print "\033[1;31mColumns not equal:"
        print t1
        print t2
        print "\033[0;0m",
        return False
    allApEq=True
    assert len(vec1)==len(vec2)
    for i in range(len(vec1)):
        if not approxSame(vec1[i],vec2[i]):
            print "\033[0;31mSignificantly different:",t1[i],"; sums:",vec1[i],",",vec2[i]
            allApEq=False
    if allApEq:
        print "No significant differences (total relative diff: "+str(approxSame.getTotalRelDiff())+"), ok."
    else:
        print "\033[1;31mSome significant differences (total relative diff: "+str(approxSame.getTotalRelDiff())+ ")!\033[0;0m"
    return allApEq

def main(fn1,fn2):
    ret=0
    opt=""
    if REL_PRECISION!=1e-6:
        opt+=" --rel-prescision="+str(REL_PRECISION)
    if ABS_PRECISION!=1e-6:
        opt+=" --abs-prescision="+str(ABS_PRECISION)
    print "\033[1;36m  compareCtsout.py"+opt+" "+fn1+" "+fn2+"\033[0;0m"
    
    # Read both files and combine into a map of key to pairs (v1, v2)
    try:
        if charEqual (fn1,fn2):
            print "ctsout.txt files are identical"
            return 0,True
        elif sumsApproxEq (fn1,fn2):
            return 0,False
        return 1,False
    # python 3000 syntax is "except IOError as e", backported to 2.6 but not always supported. Old syntax:
    except IOError, e:
        print str(e)
        return 1,False

# Test for options
def evalOptions (args):
    parser = OptionParser(usage="Usage: %prog [options] logfile1 logfile2",
            description="""Compares the simulation's ctsout.txt output files.""")
    
    parser.add_option("-R","--rel-precision",
            action="store", dest="rel_precision", type="float",
            help="Set relative precision (default: 1.0e-6)")
    parser.add_option("-A","--abs-precision",
            action="store", dest="abs_precision", type="float",
            help="Set absolute precision (default: 1.0e-6)")
    (options, others) = parser.parse_args(args=args)
    
    return options,others

if __name__ == '__main__':
    #uncomment to run unittests:
    #unittest.main()
    
    (options,others) = evalOptions (sys.argv[1:])
    if options.rel_precision:
        REL_PRECISION=options.rel_precision
    if options.abs_precision:
        ABS_PRECISION=options.abs_precision
    
    if (len(others) == 2):
        ret,ident = main (others[0],others[1])
    else:
        print "Usage: "+sys.argv[0]+" logfile1 logfile2"
        ret=-1
    sys.exit(ret)
