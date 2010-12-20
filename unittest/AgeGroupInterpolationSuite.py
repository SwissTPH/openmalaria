#!/usr/bin/python
# -*- coding: utf-8 -*-

# This file is part of OpenMalaria.
# Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licenced under the GNU General Public License version 2 or later (see file COPYING).

# This generates expected output data for the AgeGroupInterpolation unit test.
# It also plots some age-group data with interpolations.

import numpy
import matplotlib.pyplot as plt

class Interpolator(object):
    a=[]
    v=[]
    
    def __init__(self,age_lower_bounds,values):
        self.set_age_groups(age_lower_bounds,values)
    
    def set_age_groups(self,age_lower_bounds,values):
        if not len(age_lower_bounds)==len(values):
            raise AssertionError
        self.a=age_lower_bounds[:]
        self.v=values[:]
    
    def interpolate(self,t):
        i=0
        while(i+1<len(self.a) and self.a[i+1]<=t):
            i+=1
        return self.v[i]
    
    def plot(self,subplot):
        aw=[self.a[i+1]-self.a[i] for i in range(len(self.a)-1)]
        aw.append(90-self.a[len(self.a)-1])
        rects=subplot.bar(self.a, self.v, aw, 0,color='lightgreen')

class LinearInterpolator(Interpolator):
    def set_age_groups(self,age_lower_bounds,values):
        if not len(age_lower_bounds)==len(values):
            raise AssertionError
        if not age_lower_bounds[0] == 0:
            raise AssertionError
        # age group mid-points plus extra ends:
        self.a=[]
        lowb=[0]
        lowb.extend(age_lower_bounds)
        lowb.append(90) #max age — some point required
        lowb.append(1e10000) # i.e. inf
        for i in range(1,len(lowb)):
            self.a.append(0.5*(lowb[i-1]+lowb[i]))
        self.v=[values[0]]
        self.v.extend(values)
        self.v.append(values[len(values)-1])
        print "Lengths:",len(self.a),len(self.v)
    
    def interpolate(self,t):
        i=0
        while(self.a[i]<=t):
            i+=1
        if( self.a[i]-self.a[i-1] == 0.0 ):
            return self.v[i-1]
        return self.v[i-1] + (t-self.a[i-1])/(self.a[i]-self.a[i-1])*(self.v[i]-self.v[i-1])
    
    def plot(self,subplot):
        lines=subplot.plot(self.a, self.v, color='blue')


# Uses standard input data for human availability to mosquitoes:
#age group lower bounds:
ageLB=[
    0,1,2,3,4,
    5,6,7,8,9,
    10,11,12,13,14,
    15,20,20
]
#age group values:
values=[
    0.225940909648,0.286173633441,0.336898395722,0.370989854675,0.403114915112,
    0.442585112522,0.473839351511,0.512630464378,0.54487872702,0.581527755812,
    0.630257580698,0.663063362714,0.702417432755,0.734605377277,0.788908765653,
    0.839587932303,1.0,1.0
]

const=Interpolator(ageLB,values)
linear=LinearInterpolator(ageLB,values)

fig = plt.figure(1)
p = fig.add_subplot(2,1,1)
p.set_title('Human availability to mosquitoes by age')
#p.set_xlabel('age (years)')
p.set_ylabel('availability relative to an adult')

const.plot(p)
linear.plot(p)

#weight values:
values=[
    13.9856718,18.30372108,21.745749,24.25753512,26.06595444,
    28.48868784,30.84202788,33.48638244,35.20335432,37.19394024,
    40.1368962,42.00539916,44.53731348,46.77769728,49.48396092,
    54.36,60.0,60.0
]

const=Interpolator(ageLB,values)
linear=LinearInterpolator(ageLB,values)

p = fig.add_subplot(2,1,2)
p.set_title('Human mean weight by age')
p.set_xlabel('age (years)')
p.set_ylabel('mass (kg)')

const.plot(p)
linear.plot(p)

#p.set_xticks(ind+width)
#p.set_xticklabels( ('G1', 'G2', 'G3', 'G4', 'G5') )

#p.legend( (rectBars[0]), ('none') )


# Second plot: NMF data
#age group lower bounds:
ageLB=[ 0,5,10,15,60 ]
#age group values:
values=[ 6.08,3.81,2.62,4.05,5.41 ]

const=Interpolator(ageLB,values)
linear=LinearInterpolator(ageLB,values)

fig = plt.figure(2)
p = fig.add_subplot(1,1,1)
p.set_title('Interpolation methods for non-malaria fever incidence by age')
p.set_xlabel('age (years)')
p.set_ylabel('NMFs per person per year')

const.plot(p)
linear.plot(p)


# Generate sample data for unittest:
x=[15.2,18.09,7.0, 2.5, 0.0,20.0,900.0,62.0]

print "Age points:",x
print "Piecewise constant:",[const.interpolate(t) for t in x]
print "Linear interpolation:",[linear.interpolate(t) for t in x]


# Third plot: CMF data
#age group lower bounds:
ageLB=[
    0, 0.25, 0.75, 1.5,
    2.5, 3.5, 4.5, 7.5,
    12.5, 15
]
#age group values:
values=[
    0.09189, 0.0810811, 0.0648649, 0.0689189,
    0.0675676, 0.0297297, 0.0459459, 0.0945946,
    0.1243243, 0.1378378
]

const=Interpolator(ageLB,values)
linear=LinearInterpolator(ageLB,values)

fig = plt.figure(3)
p = fig.add_subplot(1,1,1)
p.set_title('Interpolation methods for case fatality rate by age')
p.set_xlabel('age (years)')
p.set_ylabel('in-patient severe malaria CFR per person per year')

const.plot(p)
linear.plot(p)


# Forth plot: CMF data
#age group lower bounds:
ageLB=[ 0.0, 5.0 ]
#age group values:
values=[ 0.0132, 0.005 ]

const=Interpolator(ageLB,values)
linear=LinearInterpolator(ageLB,values)

fig = plt.figure(4)
p = fig.add_subplot(1,1,1)
p.set_title('Interpolation methods for neurological sequelae by age')
p.set_xlabel('age (years)')
p.set_ylabel('probability of severe malaria in-patients developing sequelae')

const.plot(p)
linear.plot(p)

plt.show()
