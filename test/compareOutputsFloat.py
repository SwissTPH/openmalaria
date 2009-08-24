#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import string
import math
import unittest

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
def approx_equal (a,b, relPrecision, absPrecision):
    if isinf(a) or isinf(b):
        return False
    tolerance = relPrecision * max(math.fabs(a),math.fabs(b))
    tolerance = max(tolerance, absPrecision)
    return (math.fabs(a-b) < tolerance)

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


def main(*args):
  if (4!=len(args)):
    print "Usage: compareLogs.py logfile1 logfile2 n_prev_lines"
    return 1
  
  print "Comparing "+args[1]+" with "+args[2]
  file1=open(args[1], 'r')
  file2=open(args[2], 'r')
  pre_lines=int(args[3])+1
  line_count=0
  prev_lines1=[""]*pre_lines
  prev_lines2=[""]*pre_lines
  for line1 in file1:
    line_count+=1
    
    line_items1=string.split(line1)
    h_id1=int(line_items1[1])
    loc_id1=line_items1[2]
    value1=float(line_items1[3])
    
    line2=file2.readline()
    line_items2=string.split(line2)
    h_id2=int(line_items2[1])
    loc_id2=line_items2[2]
    value2=float(line_items2[3])
    
    if (loc_id1!=loc_id2) or (h_id1!=h_id2):
      print "Different summary outputs {0}:".format(line_count)
      print '-',line1,
      print '+',line2,
      return 2
    
    # Compare with relative precision.
    if not approx_equal_6 (value1, value2):
      print "Different values output on line {0}".format(line_count)
      print '-',line1,
      print '+',line2,
      for j in range(0, pre_lines-1):
        print '-',prev_lines1[j],
        print '+',prev_lines2[j],
      return 1
    
    if(line_count % 100000 == 0):
      print (line_count)
    for i in range(pre_lines-1, 0, -1):
      prev_lines1[i]=prev_lines1[i-1]
      prev_lines2[i]=prev_lines2[i-1]
    prev_lines1[0]=line1
    prev_lines2[0]=line2
  
  if (file2.readline() != ""):
    print "file {0} has more lines than {1}".format(args[2],args[1])
    return 3
  print "No differences, ok..."
  return 0

if __name__ == '__main__':
    #uncomment to run unittests:
    #unittest.main()
    sys.exit(main(*sys.argv))
