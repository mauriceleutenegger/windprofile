#!/usr/bin/env python2.5

import numpy as np
import matplotlib.pyplot as pl

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
u = np.arange (1000) * 0.001 * 2./3.
u = u[1:] # chop off u = 0
umax = u[-1]

# parameters
taustar = (0., 0.3, 1., 3., 10.)
q = 0.
beta = 1.
h = 0.
isNumerical = 0
isAnisotropic = 0
isRosseland = 0

# calculate fractional flux
#aatlist = []
fluxList = []
tauRlist = []
ls = ['-', '--', '-.', ':', '-']
co = ['k','k','k','k','grey']
for i in range (len (taustar)) :
    flux = windabsorption.FractionalWindEmission \
        (u, q, taustar[i], beta, h, isNumerical, isAnisotropic, isRosseland)
# u is an array (taustar q, u0, beta, h) are scalars, the rest are (bool) int
    fluxList.append (flux) 
    # renormalized each flux to make it a CD:
    flux = flux / np.max (flux)
    pl.plot (u, flux, c=co[i], ls=ls[i])
    tauR = taustar[i] * np.log (1. / (1. - u))
    tauRlist.append (tauR)
pl.xlabel ('$u$')
pl.ylabel ('$L_\lambda (u) / L_\lambda (u_0)$')
legendstring = (r'$\tau_*$ = 0', '0.3', '1', '3', '10')
leg = pl.legend (legendstring, loc='lower right', handlelength=4)
leg.draw_frame (False)
pl.xlim ([0,umax])
pl.savefig ('CD.eps')
