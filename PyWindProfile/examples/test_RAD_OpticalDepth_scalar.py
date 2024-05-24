#!/usr/bin/env python
from __future__ import print_function # for python 2 backwards compatibility

import numpy as np
import matplotlib.pyplot as pl
import PyWindProfile as wp

#p = np.arange (1.01,3.,0.01)
p = np.arange (1.01,3.,1.)
z = np.arange (-10.0, 10.001, 0.001)

#tau0 = 1.0
tau0 = 2.5e-2
EVI = 727.395
E3G = 727.049
gamma = 0.1367
deltaE = (EVI - E3G) / EVI
gammaE = gamma / EVI
vinfty_kms = 2500
speed_of_light = 2.9979e5
vinfty = vinfty_kms / speed_of_light
beta = 1.0

tscalar = np.zeros ((p.size, z.size))
for i in range (p.size) :
    for j in range (z.size) :
        tscalar[i,j] =  wp.RAD_OpticalDepth (p[i], z[j], tau0, beta,
                                             deltaE, gamma, vinfty)

#itest = 99
itest = 1

print (p[itest])

pl.semilogy (z, tscalar[itest,:], label='p={}'.format (p[itest]))
pl.xlabel ('z')
pl.ylabel ('t')
pl.legend ()
pl.savefig ('RAD_t_p2.png')
