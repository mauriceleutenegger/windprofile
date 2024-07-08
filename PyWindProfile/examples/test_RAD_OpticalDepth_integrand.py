#!/usr/bin/env python
from __future__ import print_function # for python 2 backwards compatibility

import numpy as np
import matplotlib.pyplot as pl
import PyWindProfile as wp

#p = np.arange (1.01,3.,0.01)
p = np.arange (2.0,10.1,1.)
zampl = 100.0
dzampl = 10.0
zintampl = 100.0
z0 = np.arange (-zampl,zampl+0.01, dzampl)
#z = np.arange (-10.0, 10.001, 0.001)

#tau0 = 1.0
#tau0 = 2.5e-2
#tau0 is not used for the integrand
EVI = 727.395
E3G = 727.049
gamma = 0.1367
deltaE = (EVI - E3G) / EVI
gammaE = gamma / EVI
vinfty_kms = 2500
speed_of_light = 2.9979e5
vinfty = vinfty_kms / speed_of_light
beta = 1.0

#tscalar = np.zeros ((p.size, z.size))
zlist = []
ilist = []
for i in range (p.size) :
    thiszlist = []
    thisilist = []
    for j in range (z0.size) :
        z = np.arange (z0[j], zintampl+0.001, 0.01)
        integrand = wp.RAD_OpticalDepth_integrand (z, p[i], z0[j], beta,
                                             deltaE, gammaE, vinfty)
        thiszlist.append (z)
        thisilist.append (integrand)
    zlist.append (thiszlist)
    ilist.append (thisilist)
#itest = 99
itest = -1

print (p[itest])

for j in range (z0.size) :
    z = zlist[itest][j]
    integrand = ilist[itest][j]
    pl.semilogy (z, integrand, label='p={}, z0={}'.format (p[itest], z0[j]))
pl.xlabel ('z')
pl.ylabel ('t')
pl.legend ()
pl.savefig ('RAD_i_p10.png')
