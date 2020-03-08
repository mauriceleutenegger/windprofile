#!/usr/bin/env python
from __future__ import print_function # for python 2 backwards compatibility

import numpy as np
import matplotlib.pyplot as pl
import PyWindProfile as wp

z = np.arange (1.01, 10.001, 0.01)
p = np.array ([0., 1.0, 2., 3.])

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

for i in range (p.size) :
    labelstring = 'p = {}'.format (p[i])
    pl.plot (z, t[i,:], label=labelstring)
    pl.ylim (0., 2.)
pl.xlabel ('z ($R_*$)')
pl.ylabel ('t (p,z)')
pl.legend ()
pl.savefig ('t_z.eps')


