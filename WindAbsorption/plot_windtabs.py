#!/usr/bin/env python2.5
"""
    Plot output of generate_windtabs.py and porous_windtabs.py
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
from pyfits import *

files = ('tau_transmission.fits', 'tau_transmission1.fits',\
             'tau_transmission3.fits', 'tau_transmission10.fits')
legendtext = ('h = 0', '1', '3', '10')

for file in files :
    data = getdata (file)
    tau = data.field ('tau')
    transmission = data.field ('transmission')
    semilogy (tau, transmission)

xlabel (r'$\tau_*$')
ylabel ('Transmission')
#xlim (0,40)
xlim (0,1)
ylim (5.e-3,1)
legend (legendtext)
savefig ('porous_windtabs.eps')
