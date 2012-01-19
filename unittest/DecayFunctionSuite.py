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

# This isn't actually used to generate unit-test values, but just to graph
# decay functions.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import cnames
import string
import math as m


fig = plt.figure(1)

p = fig.add_subplot(2,1,1)
p.set_title("Decay functions")
p.set_xlabel('time (years)')
p.set_ylabel('decayed value (L=10)')

L = 10
log2=m.log(2.0)
x=np.arange(0,30,0.025)

lines=list()
labels=list()

y=[m.exp(-t/L * log2) for t in x]
pl=p.plot(x,y,color='cyan')
lines.append(pl)
labels.append('exponential')

k=0.5

y=[m.exp(-(t/L)**k * log2) for t in x]
pl=p.plot(x,y,color='green')
lines.append(pl)
labels.append('weibull (k=0.5)')

y=[1./(1 + (t/L)**k) for t in x]
pl=p.plot(x,y,color='red')
lines.append(pl)
labels.append('hill (k=0.5)')

k=2

y=[m.exp(-(t/L)**k * log2) for t in x]
pl=p.plot(x,y,color='violet')
lines.append(pl)
labels.append('weibull (k=2)')

y=[1./(1 + (t/L)**k) for t in x]
pl=p.plot(x,y,color='orange')
lines.append(pl)
labels.append('hill (k=2)')

plt.legend(lines,labels,'upper right')

p = fig.add_subplot(2,1,2)
p.set_xlabel('time (years)')
p.set_ylabel('decayed value (L=20)')

L = 20
k=0.5

lines=list()
labels=list()

y=[(1-t/L if t<L else 0) for t in x]
pl=p.plot(x,y,color='brown')
lines.append(pl)
labels.append('linear')

y=[m.exp( k - k / (1 - (t/L)**2) ) if t<L else 0 for t in x]
pl=p.plot(x,y,color='green')
lines.append(pl)
labels.append('smooth-compact (k=0.5)')

k=1
y=[m.exp( k - k / (1 - (t/L)**2) ) if t<L else 0 for t in x]
pl=p.plot(x,y,color='blue')
lines.append(pl)
labels.append('smooth-compact (k=1)')

k=2
y=[m.exp( k - k / (1 - (t/L)**2) ) if t<L else 0 for t in x]
pl=p.plot(x,y,color='lightblue')
lines.append(pl)
labels.append('smooth-compact (k=2)')

plt.legend(lines,labels,'upper right')

plt.show()
