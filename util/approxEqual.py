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

import math
import unittest

# Could use math.isinf, but it's not present in python 2.5
def isinf (a):
    return (a-a) != 0
def isnan (a):
    return a != a

# Object to test whether two values are equal to a certain precision
class ApproxEqual (object):
    def __init__(self, relPrec, absPrec):
        self.relPrecision = relPrec
        self.absPrecision = absPrec
        self.sumRelDiffs = 0.0
    
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
    def __call__ (self,a,b):
        if isinf(a) or isinf(b):
            relDiff = 1e10000 * 0.0 # NaN
        else:
            tolerance = self.relPrecision * max(math.fabs(a),math.fabs(b))
            tolerance = max(tolerance, self.absPrecision)
            relDiff = math.fabs(a-b) / tolerance
        self.sumRelDiffs += relDiff
        return (relDiff < 1.0)
    
    def getTotalRelDiff(self):
        # this value is for example increased (from 0) by 1/2 by every
        # difference by an order of 2, as long as the values are larger
        # than absPrecision. Each difference should add less than 1.
        return self.sumRelDiffs * self.relPrecision

# As ApproxEqual, but considers two NaNs or two pos./neg. infs "the same"
class ApproxSame (ApproxEqual):
    def __init__(self, relPrec, absPrec):
        ApproxEqual.__init__(self, relPrec, absPrec)
    
    def __call__ (self,a,b):
        if isnan(a) and isnan(b):
            return True
        elif isinf(a) and isinf(b):
            if copysign(1.0,a)==copysign(1.0,b):
                return True
            else:
                self.sumRelDiffs += 1e10000
                return False
        else:
            return ApproxEqual.__call__(self,a,b)

class TestApproxEqual (unittest.TestCase):
    approxEqual6 = ApproxEqual(1e-6, 1e-6)
    
    def testNaN (self):
        Inf = 1e10000
        NaN = Inf * 0.0
        self.assertTrue (not self.approxEqual6 (NaN, 1))
        self.assertTrue (not self.approxEqual6 (NaN, 0))
        self.assertTrue (not self.approxEqual6 (NaN, NaN))
        self.assertTrue (not self.approxEqual6 (NaN, Inf))
    
    def testInf (self):
        Inf = 1e10000
        self.assertTrue (not self.approxEqual6 (Inf, 1))
        self.assertTrue (not self.approxEqual6 (Inf, 0))
        self.assertTrue (not self.approxEqual6 (Inf, Inf))
        self.assertTrue (not self.approxEqual6 (Inf, -Inf))
    
    def testZero (self):
        self.assertTrue (not self.approxEqual6 (0, 1e-6))
        self.assertTrue (    self.approxEqual6 (0, 1e-7))
    
    def testRegular (self):
        self.assertTrue (not self.approxEqual6 (1, 0))
        self.assertTrue (not self.approxEqual6 (1, 0.999999))
        self.assertTrue (    self.approxEqual6 (1, 0.9999995))
        self.assertTrue (not self.approxEqual6 (1000000, 999999))
        self.assertTrue (    self.approxEqual6 (1000000, 999999.5))
        # these are considered equal because of absolute precision limitation rather than relative:
        self.assertTrue (    self.approxEqual6 (0.000001, 0.0000005))
        # this is roughly on the verge of what isn't considered equal:
        self.assertTrue (not self.approxEqual6 (0.000001, 0.000002))
        
        # if we only want to test relative precision:
        approxEqualRel = ApproxEqual(1e-6, 0)
        self.assertTrue (not approxEqualRel (0.000001, 0.000000999999))
        self.assertTrue (    approxEqualRel (0.000001, 0.0000009999995))
        approxEqual = ApproxEqual(0, 1e-11)
        self.assertTrue (approxEqual (approxEqualRel.sumRelDiffs, 1.49999999994))


if __name__ == '__main__':
    # run unittests:
    unittest.main()
