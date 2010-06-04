#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import string
import math
import unittest
from optparse import OptionParser

totalRelDiff = 0.0

REL_PRECISION=1e-6
ABS_PRECISION=1e-6

# Could use math.isinf, but it's not present in python 2.5
def isinf (a):
    return (a-a) != 0

# Careful with NaN, +/- inf and 0 values! Note: inf == inf
# Check a and b are approximately equal. Return true if:
#   a equals b to at least log10(relPrecision) significant figures
#   or at least log10(absPrecision) decimal places.
# This should work for small and large values, when one is zero, and when
# either is infinite or an NaN.
# 
# What it is limited by, is testing of small values with relative precision
# when it should also allow for large values being rounded to zero.
# Possibly considering the case where a or b is zero explicitly would help, but
# might lead to confusing outcomes.
def approx_equal (a,b, relPrecision=1e-6, absPrecision=1e-6):
    if isinf(a) or isinf(b):
        relDiff = 1e1000-1e1000 # NaN
    else:
        tolerance = relPrecision * max(math.fabs(a),math.fabs(b))
        tolerance = max(tolerance, absPrecision)
        relDiff = math.fabs(a-b) / tolerance
    global totalRelDiff
    totalRelDiff += relDiff
    return (relDiff < 1.0)

def approx_equal_6 (a,b):
    return approx_equal (a, b, 1e-6, 1e-6)

class TestApproxEqual (unittest.TestCase):
    def testNaN (self):
        Inf = 1e10000
        NaN = Inf * 0
        self.assert_ (not approx_equal_6 (NaN, 1))
        self.assert_ (not approx_equal_6 (NaN, 0))
        self.assert_ (not approx_equal_6 (NaN, NaN))
        self.assert_ (not approx_equal_6 (NaN, Inf))
    
    def testInf (self):
        Inf = 1e10000
        self.assert_ (not approx_equal_6 (Inf, 1))
        self.assert_ (not approx_equal_6 (Inf, 0))
        self.assert_ (not approx_equal_6 (Inf, Inf))
        self.assert_ (not approx_equal_6 (Inf, -Inf))
    
    def testZero (self):
        self.assert_ (not approx_equal_6 (0, 1e-6))
        self.assert_ (    approx_equal_6 (0, 1e-7))
    
    def testRegular (self):
        self.assert_ (not approx_equal_6 (1, 0))
        self.assert_ (not approx_equal_6 (1, 0.999999))
        self.assert_ (    approx_equal_6 (1, 0.9999995))
        self.assert_ (not approx_equal_6 (1000000, 999999))
        self.assert_ (    approx_equal_6 (1000000, 999999.5))
        # these are considered equal because of absolute precision limitation rather than relative:
        self.assert_ (    approx_equal_6 (0.000001, 0.0000005))
        # this is roughly on the verge of what isn't considered equal:
        self.assert_ (not approx_equal_6 (0.000001, 0.000002))
        # if we only want to test relative precision:
        self.assert_ (not approx_equal (0.000001, 0.000000999999,  1e-6, 0))
        self.assert_ (    approx_equal (0.000001, 0.0000009999995, 1e-6, 0))


class ValIdentifier:
    def __init__(self,s,g,m):
        self.survey = s
        self.group = g
        self.measure = m
    def __eq__(self,other):
        return (self.survey == other.survey) and (self.group == other.group) and (self.measure == other.measure)
    def __hash__(self):
        return self.survey.__hash__() ^ self.group.__hash__() ^ self.measure.__hash__()

class TestValIdentifier (unittest.TestCase):
    def setUp(self):
        self.a1 = ValIdentifier(2,4,0);
        self.a2 = ValIdentifier(2,4,0);
        self.a3 = ValIdentifier(2,2,0);
        self.b1 = ValIdentifier(0,"abc",5);
        self.b2 = ValIdentifier(0,"abc",5);
    def testEq (self):
        self.assert_ (self.a1.survey == self.a2.survey)
        self.assert_ (self.a1 == self.a2)
        self.assert_ (self.b1 == self.b2)
        self.assert_ (self.a1 != self.a3)
    def testHash (self):
        self.assert_ (self.a1.__hash__ == self.a2.__hash__)
        self.assert_ (self.a1.__hash__ != self.a3.__hash__) # actually, could sometimes be
        self.assert_ (self.b1.__hash__ == self.b2.__hash__)


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

def ReadEntries (fname):
    values=dict()
    fileObj = open(fname, 'r')
    for line in fileObj:
        items=string.split(line)
        if (len(items) != 4):
            print "expected 4 items on line; found (following line):"
            print line
            continue
            
        key=ValIdentifier(int(items[0]),items[1],int(items[2]))
        values[key]=float(items[3])
    return values

def main(fn1,fn2,maxDiffsToPrint=6):
    """Takes names of the two files to compare and optionally an argument describing
the maximum number of differences to print directly (note: order is not intuitive).
Returns a tuple ret,ident; ret is 0 if test passes (output considered near-enough equal),
ident is 1 if files are binary-equal."""
    ret=0
    opt=""
    if REL_PRECISION!=1e-6:
        opt+=" --rel-prescision="+str(REL_PRECISION)
    if ABS_PRECISION!=1e-6:
        opt+=" --abs-prescision="+str(ABS_PRECISION)
    print "compareOutputsFloat.py"+opt+" "+fn1+" "+fn2+" "+str(maxDiffsToPrint)
    
    # Read both files and combine into a map of key to pairs (v1, v2)
    try:
        if charEqual (fn1,fn2):
            print "Files are identical"
            return 0,True
        print "Files aren't binary-equal"
        
        values1=ReadEntries(fn1)
        values2=ReadEntries(fn2)
    # python 3000 syntax is "except IOError as e", backported to 2.6 but not always supported. Old syntax:
    except IOError, e:
        print str(e)
        return 1,False
    values=dict()
    for (k,v1) in values1.iteritems():
        v2=None
        if (k in values2):
            v2=values2[k]
            del values2[k]
        values[k] = (v1,v2)
    for (k,v2) in values2.iteritems():
        values[k] = (None,v2)
    
    # Go through all values:
    numPrinted=0
    numDiffs=0
    numMissing1=0
    numMissing2=0
    perMeasureNum = dict()
    perMeasureTotal = dict()
    perMeasureNumDiff = dict()
    perMeasureDiffSum = dict()
    perMeasureDiffAbsSum = dict()
    for (k,(v1,v2)) in values.iteritems():
        if v1==None:
            numMissing1 += 1
        elif v2==None:
            numMissing2 += 1
        else:
            perMeasureNum[k.measure] = perMeasureNum.get(k.measure, 0) + 1
            perMeasureTotal[k.measure] = perMeasureTotal.get(k.measure, 0.0) + v1
            
            # Compare with relative precision
            if approx_equal (v1, v2, REL_PRECISION, ABS_PRECISION):
                continue
            
            numDiffs += 1
            # Sum up total difference per measure
            perMeasureDiffSum[k.measure]    = perMeasureDiffSum.get(k.measure,0.0)    + v2 - v1
            perMeasureDiffAbsSum[k.measure] = perMeasureDiffAbsSum.get(k.measure,0.0) + math.fabs(v2-v1)
        
        numPrinted += 1
        perMeasureNumDiff[k.measure] = perMeasureNumDiff.get(k.measure,0) + 1;
        if (numPrinted <= maxDiffsToPrint):
            print "survey "+str(k.survey)+", group "+str(k.group)+", measure "+str(k.measure)+": "+str(v1)+" -> "+str(v2)
            if (numPrinted == maxDiffsToPrint):
                print "[won't print any more line-by-line diffs]"
    
    if (numMissing1 > 0) or (numMissing2 > 0):
        print str(numMissing1) + " entries missing from first file, " + str(numMissing2) +" from second"
        ret = 3
    
    for (k.measure,absDiff) in perMeasureDiffAbsSum.iteritems():
        if not (absDiff <= 1e-6):   # handle NANs
            # standard division throws on divide-by-zero, which I don't want
            def div(x,y):
                try:
                    return x/y
                except ZeroDivisionError:
                    return 1e400 * 0 # nan
            
            diff=perMeasureDiffSum[k.measure]
            sum1=perMeasureTotal[k.measure]
            print "for measure "+str(k.measure)+": sum(1st file):"+str(sum1)+" diff/sum: "+str(div(diff,sum1))+" (abs diff)/sum: "+str(div(absDiff,sum1))+" diff/(abs diff): "+str(div(diff,absDiff))+"  num diffs/total: "+str(perMeasureNumDiff.get(k.measure,0)/perMeasureNum.get(k.measure,0))
    
    # We print total relative diff here: 1.0 should mean roughly, one parameter is twice what it should be.
    if numDiffs == 0:
        print "No significant differences (total relative diff: "+str(totalRelDiff/1.e6)+"), ok..."
        return ret,False
    else:
        print str(numDiffs)+" significant differences (total relative diff: "+str(totalRelDiff/1.e6)+ ")!"
        return 1,False

# Test for options
def evalOptions (args):
    parser = OptionParser(usage="Usage: %prog [options] logfile1 logfile2 [max different lines to print]",
            description="""Scenarios to be run must be of the form scenarioXX.xml; if any are passed on the command line, XX is substituted for each given; if not then all files of the form scenario*.xml are run as test scenarios.
You can pass options to openMalaria by first specifying -- (to end options passed from the script); for example: %prog 5 -- --print-model""")
    
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
    
    if (len(others) == 3):
        ret,ident = main (others[0],others[1],int(others[2]))
    elif (len(others) == 2):
        ret,ident = main (others[0],others[1])
    else:
        print "Usage: "+sys.argv[0]+" logfile1 logfile2 [max different lines to print]"
        ret=-1
    sys.exit(ret)
