#!/usr/bin/env python

import numpy as np
import pyfits as pf

# TauStar grid

# linear spacing for 0 < taustar < 1; 
# followed by logarithmic spacing for 1 < taustar < 1000

size01 = 1000
sizeE0E3 = 3000 
taustar01 = np.arange (size01) / float (size01)
logtaustarE0E3 = 3. * np.arange (sizeE0E3) / float (sizeE0E3)
taustarE0E3 = 10.**logtaustarE0E3
taustar = np.concatenate ([taustar01, taustarE0E3])
sizeT = size01 + sizeE0E3

#sizeK = 1000
sizeK = 50 # add in successive pieces as the calculation runs
kappaRatio = np.arange (sizeK) / 50. # run from zero to twenty in steps of .02

transmission = np.zeros ((sizeT, sizeK))

for i in range (sizeK):
    k = kappaRatio[i]
    infile = 'tau_transmission_HeII_kappa_' + str (k) + '.txt'
    t_k = np.loadtxt (infile)
    transmission[:,i] = t_k

# write to a fits file

tauCol = pf.Column (name='tau', format = 'E', array = taustar)
transmissionCol = pf.Column (name='transmission', format='E', array=transmission, dim='2')
kappaCol = pf.Column (name='kappaRatio', format='E', array = kappaRatio)

columns = pf.ColDefs ([tauCol, ])
columns2 = pf.ColDefs ([kappaCol,])
columns3 = pf.ColDefs ([transmissionCol, ])

primaryHDU = pf.PrimaryHDU (data=transmission)
tableHDU = pf.new_table (columns)
tableHDU2 = pf.new_table (columns2)

outfile = 'tau_transmission_HeII.fits'
hdulist = pf.HDUList ([primaryHDU, tableHDU, tableHDU2])
hdulist.writeto (outfile, clobber=True)
