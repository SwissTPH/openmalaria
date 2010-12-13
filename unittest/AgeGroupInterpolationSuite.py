#!/usr/bin/python
# -*- coding: utf-8 -*-

# This file is part of OpenMalaria.
# Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licenced under the GNU General Public License version 2 or later (see file COPYING).

# This generates expected output data for the AgeGroupInterpolation unit test.

# Uses standard input data for human availability to mosquitoes:
#age group lower bounds:
a=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,20]
#age group values:
v=[0.225940909648,0.286173633441,0.336898395722,0.370989854675,0.403114915112,0.442585112522,0.473839351511,0.512630464378,0.54487872702,0.581527755812,0.630257580698,0.663063362714,0.702417432755,0.734605377277,0.788908765653,0.839587932303,1.0,1.0]

def piecewiseConst(t):
    i=0
    while(i+1<len(a) and a[i+1]<=t):
        i+=1
    return v[i]

# age group mid-points plus extra ends:
b=[]
l=0
for x in a:
    b.append(0.5*(l+x))
    l=x
b.append(1e10000) # i.e. inf
v.append(v[len(v)-1]) # repeat last value

def lerp(t):
    i=0
    while(b[i]<=t):
        i+=1
    if( b[i]-b[i-1] == 0.0 ):
        return v[i-1]
    return v[i-1] + (t-b[i-1])/(b[i]-b[i-1])*(v[i]-v[i-1])

x=[15.2,18.09,7.0, 2.5, 0.0,20.0,900.0]

print "Age points:",x
print "Piecewise constant:",[piecewiseConst(t) for t in x]
print "Linear interpolation:",[lerp(t) for t in x]
