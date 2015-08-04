#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl

import os
HomePath = os.path.expanduser ("~")
LibPath = os.path.join (HomePath, 'lib/python/PyWindProfile/')
import sys
sys.path.append (LibPath)
import PyWindProfile as wp

#psize = 1000
zsize = 1000
#p = np.arange (psize) / 100. + 1.1 
z = np.arange (zsize) / 10. + 1.01
psize = 3
p = np.array ([0., 1.1, 2., 3.])
#z = np.array ([0., 1., 1.1, 2., 5., 10.])

t = wp.OpticalDepth2d (p,z,1.,0.,1.,1,0,0,0,0)
tHeII = wp.OpticalDepth2d (p,z,1.,0.,1.,1,0,0,0,1)
ratio = t[0,:] / tHeII[0,:]

pl.semilogx (z, t[0,:])
pl.semilogx (z, tHeII[0,:])
pl.ylim (0., 2.)
pl.savefig ('t_p0.eps')

#print "finished 2d version"

#t2 = t
#for i in range (psize) :
#    for j in range (zsize) :
#      t[i,j] = wp.OpticalDepth (p[i],z[i],1.,0.,1.,0,0,0,0)

#print "finished scalar version"
#t = wp.OpticalDepth2d (p,z,1.,0.,1.,0,0,0,0)

#print "finished 2d version again"
