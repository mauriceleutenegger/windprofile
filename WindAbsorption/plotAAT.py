#!/usr/bin/env python
"""
Docstring here...
"""
import numpy as np
import matplotlib.pyplot as pl

# font stuff
#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)


# from scipy cookbook for publication quality figures:

fig_width_pt = 245.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inches
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height =fig_width*0.75       # height in inches
fig_size = [fig_width,fig_height]
params = {'backend' : 'eps',
          'axes.labelsize' : 10,
          'font.size' : 10,
          'font.family' : 'sans-serif',
          'font.family.sans-serif' : ['Helvetica'],
          'text.fontsize' : 10,
          'legend.fontsize': 6,
          'xtick.labelsize' : 8,
          'ytick.labelsize' : 8,
          'text.usetex' : True,
          'figure.figsize' : fig_size}
pl.rcParams.update (params)

# make the (first) plot
pl.figure(1)
pl.clf()
pl.axes([0.125,0.2,0.95-0.125,0.95-0.2])

import os
HomePath = os.path.expanduser ("~")
WindAbsorptionPath = os.path.join (HomePath, 'lib/python/windabsorption')
import sys
sys.path.append (WindAbsorptionPath)
import windabsorption

# inverse radius array
u = np.arange (900) * 0.001 + 0.001
#u = np.arange (1000) * 0.001 
#u = u[1:-10]

# parameters
taustar = (0., 0.3, 1., 3., 10.)
beta = 1.
h = 0.
isNumerical = 0
isAnisotropic = 0
isRosseland = 0
kappaRatio = 20.
isHeII = 1

# calculate aat
aatlist = []
tauRlist = []
ls = ['-', '--', '-.', ':', '-']
co = ['k','k','k','k','grey']
for i in range (len (taustar)) :
    aat = windabsorption.AngleAveragedTransmission \
        (u, taustar[i], kappaRatio, beta, h, isNumerical, isAnisotropic, isRosseland, isHeII)
# u is an array (taustar q, u0, beta, h) are scalars, the rest are (bool) int
    aatlist.append (aat) 
    pl.plot (u, aat, c=co[i], ls=ls[i])
    tauR = taustar[i] * np.log (1. / (1. - u))
    tauRlist.append (tauR)
pl.xlabel ('$u$')
pl.ylabel ('$\overline{T}$')
#pl.ylabel ('T')
legendstring = (r'$\tau_*$ = 0', '0.3', '1', '3', '10')
leg = pl.legend (legendstring, loc='center right', handlelength=4)
#leg.draw_frame (False)
#pl.title ('Angle Averaged Transmission')
pl.savefig ('aat.eps')

# need to add a comparison to e^-tau
# x-axis should be radial optical depth

for i in range (len (taustar)) :
    tauR = tauRlist[i]
    aat = aatlist[i]
    pl.semilogy (tauR, aat, hold=bool(i))
pl.semilogy (tauR, np.exp (-1. * tauR), color='k', ls='--')
pl.semilogy (tauR, np.exp (-2. * tauR), color='b', ls='--')
legendstring = (r'$\tau_* = 0$', '0.3', '1', '3', '10', r'$e^{-\tau}$', 
                r'$e^{-2\tau}$')
pl.legend (legendstring, loc='upper right')
pl.xlabel (r'$\tau_R$')
pl.ylabel ('T')
pl.title ('Angle Averaged Transmission')
pl.xlim (0,4)
pl.ylim (1.e-2,1.)
pl.savefig ('aat2.eps')
