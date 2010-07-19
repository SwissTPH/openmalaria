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

import math
import unittest

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
def approxEqual (a,b, relPrecision=1e-6, absPrecision=1e-6):
    if isinf(a) or isinf(b):
        relDiff = 1e1000-1e1000 # NaN
    else:
        tolerance = relPrecision * max(math.fabs(a),math.fabs(b))
        tolerance = max(tolerance, absPrecision)
        relDiff = math.fabs(a-b) / tolerance
    global totalRelDiff
    totalRelDiff += relDiff
    return (relDiff < 1.0)

def approxEqual6 (a,b):
    return approxEqual (a, b, 1e-6, 1e-6)

class TestApproxEqual (unittest.TestCase):
    def testNaN (self):
        Inf = 1e10000
        NaN = Inf * 0
        self.assert_ (not approxEqual6 (NaN, 1))
        self.assert_ (not approxEqual6 (NaN, 0))
        self.assert_ (not approxEqual6 (NaN, NaN))
        self.assert_ (not approxEqual6 (NaN, Inf))
    
    def testInf (self):
        Inf = 1e10000
        self.assert_ (not approxEqual6 (Inf, 1))
        self.assert_ (not approxEqual6 (Inf, 0))
        self.assert_ (not approxEqual6 (Inf, Inf))
        self.assert_ (not approxEqual6 (Inf, -Inf))
    
    def testZero (self):
        self.assert_ (not approxEqual6 (0, 1e-6))
        self.assert_ (    approxEqual6 (0, 1e-7))
    
    def testRegular (self):
        self.assert_ (not approxEqual6 (1, 0))
        self.assert_ (not approxEqual6 (1, 0.999999))
        self.assert_ (    approxEqual6 (1, 0.9999995))
        self.assert_ (not approxEqual6 (1000000, 999999))
        self.assert_ (    approxEqual6 (1000000, 999999.5))
        # these are considered equal because of absolute precision limitation rather than relative:
        self.assert_ (    approxEqual6 (0.000001, 0.0000005))
        # this is roughly on the verge of what isn't considered equal:
        self.assert_ (not approxEqual6 (0.000001, 0.000002))
        # if we only want to test relative precision:
        self.assert_ (not approxEqual (0.000001, 0.000000999999,  1e-6, 0))
        self.assert_ (    approxEqual (0.000001, 0.0000009999995, 1e-6, 0))


if __name__ == '__main__':
    # run unittests:
    unittest.main()
