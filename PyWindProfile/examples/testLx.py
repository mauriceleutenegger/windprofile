#!/usr/bin/env python
from __future__ import print_function # for python 2 backwards compatibility

import numpy as np
import matplotlib.pyplot as pl
import PyWindProfile as wp

# set xgrid
dx = 1.e-3
x = np.arange (-1., 1.+0.1*dx, dx)

# set parameters
q = 0.
U0 = 0.65
Umin0 = 0.
beta = 1.
TauStar = 2.
h = 0.
Tau0Star = 0.
betaSobolev = 0.
kappaRatio = 0.
isOpticallyThick = 0
isNumerical = 0
isAnisotropic = 0
isRosseland = 0
isExpansion = 0

# calculate line profile
flux = wp.Lx\
    (x, q, U0, Umin0, beta, TauStar, h, Tau0Star, betaSobolev, kappaRatio, \
     isOpticallyThick, isNumerical, isAnisotropic, isRosseland, \
     isExpansion, 0);

# plot
pl.plot (x,flux,c='k')
pl.xlabel ('x')
pl.ylabel ('flux')
titlestr = '$\\tau_* = {}, u_0 = {}$'.format (TauStar, U0)
pl.title (titlestr)
pl.savefig ('lx.eps')

