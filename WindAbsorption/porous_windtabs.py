#!/usr/bin/env python2.5
"""
    Generate a FITS file giving T(tau) for porous winds with certain parameters.
"""

"""
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA */
"""
from numpy import *
from pylab import *
import pyfits
import os
HomePath = os.path.expanduser ("~")
WindAbsorptionPath = os.path.join (HomePath, 'lib/python/windabsorption')
import sys
sys.path.append (WindAbsorptionPath)
import windabsorption
#WindAbsorption

# just make it a little easier to write out the data
def writeFits (taustar, transmission, filename) :
    tauCol = pyfits.Column (name='tau', format = 'E', array = taustar)
    transmissionCol = pyfits.Column (name='transmission', format='E', 
                                     array=transmission)
    columns = pyfits.ColDefs ([tauCol, transmissionCol])
    tableHDU = pyfits.new_table (columns)
    tableHDU.writeto (filename)



# linear spacing for 0 < taustar < 1; 
# followed by logarithmic spacing for 1 < taustar < 1000
size01 = 1000
sizeE0E3 = 3000
taustar01 = arange (size01) / float (size01)
logtaustarE0E3 = 3. * arange (sizeE0E3) / float (sizeE0E3)
taustarE0E3 = 10.**logtaustarE0E3
taustar = concatenate ([taustar01, taustarE0E3])

# parameters
q = 0.
u0 = 2./3.
beta = 1.
h = 1.
isNumerical = 0
isAnisotropic = 0
isRosseland = 0

# calculate transmission
transmission1 = windabsorption.WindAbsorption (taustar, q, u0, beta, h, isNumerical, isAnisotropic, isRosseland)
# taustar is an array, (q, u0, beta, h) are scalars, the rest are (bool) int

# write to a fits file
writeFits (taustar, transmission1, 'tau_transmission1.fits')

h = 3.
transmission3 = windabsorption.WindAbsorption (taustar, q, u0, beta, h, isNumerical, isAnisotropic, isRosseland)
writeFits (taustar, transmission3, 'tau_transmission3.fits')

h = 10.
transmission10 = windabsorption.WindAbsorption\
    (taustar, q, u0, beta, h, isNumerical, isAnisotropic, isRosseland)
writeFits (taustar, transmission10, 'tau_transmission10.fits')


