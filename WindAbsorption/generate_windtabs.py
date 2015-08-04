#!/usr/bin/env python
"""
    Generates a FITS file giving T(tau) for a non-porous wind with certain parameters.
    This file is then used by the XSPEC local model windtabs.
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

import numpy as np
import pyfits as pf
import os
# Because python only looks in ~/lib/python
HomePath = os.path.expanduser ("~")
WindAbsorptionPath = os.path.join (HomePath, 'lib/python/windabsorption')
import sys
sys.path.append (WindAbsorptionPath)
import windabsorption
#WindAbsorption

# linear spacing for 0 < taustar < 1; 
# followed by logarithmic spacing for 1 < taustar < 1000
size01 = 1000
sizeE0E3 = 3000
taustar01 = np.arange (size01) / float (size01)
logtaustarE0E3 = 3. * np.arange (sizeE0E3) / float (sizeE0E3)
taustarE0E3 = 10.**logtaustarE0E3
taustar = np.concatenate ([taustar01, taustarE0E3])

# parameters
q = 0.
u0 = 2./3.
beta = 1.
h = 0.
isNumerical = 0
isAnisotropic = 0
isRosseland = 0

# calculate transmission
transmission = windabsorption.WindAbsorption (taustar, q, u0, beta, h, isNumerical, isAnisotropic, isRosseland)
# taustar is an array, (q, u0, beta, h) are scalars, the rest are (bool) int

# write to a fits file

tauCol = pf.Column (name='tau', format = 'E', array = taustar)
transmissionCol = pf.Column (name='transmission', format='E', array=transmission)

columns = pf.ColDefs ([tauCol, transmissionCol])

tableHDU = pf.new_table (columns)

tableHDU.writeto ('tau_transmission.fits')
