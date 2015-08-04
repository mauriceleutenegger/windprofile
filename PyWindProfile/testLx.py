#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl
import os
HomePath = os.path.expanduser ("~")
LibPath = os.path.join (HomePath, 'lib/python/PyWindProfile/')
import sys
sys.path.append (LibPath)
import PyWindProfile as wp

# set xgrid
xsize=1000.
x = np.arange (2*xsize+1) / xsize - 1. 

# set parameters
q = 0.
U0 = 0.65
Umin0 = 0.
#Umin02 = 0.2
beta = 1.
TauStar = 2.
h = 0.
Tau0Star = 0.
betaSobolev = 0.
isOpticallyThick = 0
isNumerical = 0
isAnisotropic = 0
isRosseland = 0
isExpansion = 0

#kappaRatio = 3.

# calculate line profile
flux = wp.Lx\
    (x, q, U0, Umin0, beta, TauStar, h, Tau0Star, betaSobolev, 0., \
     isOpticallyThick, isNumerical, isAnisotropic, isRosseland, \
     isExpansion, 0);
flux_kr1 = wp.Lx\
    (x, q, U0, Umin0, beta, TauStar, h, Tau0Star, betaSobolev, 1.,\
     isOpticallyThick, isNumerical, isAnisotropic, isRosseland, \
     isExpansion, 1);
flux_kr3 = wp.Lx\
    (x, q, U0, Umin0, beta, TauStar, h, Tau0Star, betaSobolev, 3.,\
     isOpticallyThick, isNumerical, isAnisotropic, isRosseland, \
     isExpansion, 1);

#flux5 = wp.Lx\
#    (x, q, U0, Umin02, beta, TauStar, h, Tau0Star, betaSobolev,\
#         isOpticallyThick, isNumerical, isAnisotropic, isRosseland,\
#         isExpansion);
# renormalize
print flux.sum ()
print flux_kr1.sum ()
print flux_kr3.sum ()
flux /= flux.max ()
flux_kr1 /= flux_kr1.max ()
flux_kr3 /= flux_kr3.max ()

# plot
pl.plot (x,flux,c='k', label = 'no He+')
pl.plot (x,flux_kr1, c='r', label = 'kappaRatio = 1')
pl.plot (x,flux_kr3, c='b', label = 'kappaRatio = 3')
pl.xlabel ('x')
pl.ylabel ('flux')
pl.legend (loc='upper right')
pl.title ('$\\tau_* = '+str(TauStar)+', u_0 = 0.65 R_*$')
pl.savefig ('lx.eps')

#pl.clf ()
# ratio plot
#ratio = flux_kr1 / flux
#pl.plot (x, ratio, c='k')
#pl.xlabel ('x')
#pl.ylabel ('ratio')
#pl.title ('$\\tau_* = 3, u_0 = 0.65 R_*$')
#pl.savefig ('lx_ratio.eps')


#np.savetxt ('lx.txt', np.transpose ([x, flux, flux5]), fmt='%g')
