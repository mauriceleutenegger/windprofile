#!/usr/bin/env python
from __future__ import print_function # for python 2 backwards compatibility

import numpy as np
import matplotlib.pyplot as pl
import PyWindProfile as wp

# calculate t(p,z) for a large grid of points that shows
# how much faster it is to calculate using the 2d-array returning function
# vs. the scalar function in a loop

p = np.arange (0.,3.,0.01)
z = np.arange (1.01, 10.001, 0.001)

taustar = 1.0
h = 0.0
beta = 1.0
numerical = 1
anisotropic = 0
rosseland = 0
expansion = 0
HeII = 0

t = wp.OpticalDepth2d (p, z, taustar, h, beta, numerical, anisotropic,
                       rosseland, expansion, HeII)
print ("finished 2d version")
tscalar = t
psize = p.size
zsize = z.size
for i in range (psize) :
    for j in range (zsize) :
      tscalar[i,j] = wp.OpticalDepth (p[i], z[i], taustar, h, beta,
                                      numerical, anisotropic, rosseland,
                                      expansion)
print ("finished scalar version")

t = wp.OpticalDepth2d (p, z, taustar, h, beta, numerical, anisotropic,
                       rosseland, expansion, HeII)
print ("finished 2d version again")

