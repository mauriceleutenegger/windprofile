/***************************************************************************
    RAD_OpticalDepthZ.h   - Calculates the optical depth to X-rays from 
                           Resonant Auger Destruction (RAD) for a stellar
                           wind along a ray p from point z.

                             -------------------
    begin				: May 2024
    copyright			: (C) 2024 by Gabe Grell and Maurice Leutenegger
    email				: gabriel.j.grell@nasa.gov maurice.a.leutenegger@nasa.gov
 ***************************************************************************/
 /* This program is free software; you can redistribute it and/or modify
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

#include "RAD_OpticalDepth.h"

RAD_OpticalDepthZ::RAD_OpticalDepthZ (RAD_OpticalDepth* A)
  : Integral (), itsRADOD (A)
{
  return;
}

// z1 < z2
Real RAD_OpticalDepthZ::getOpticalDepth (Real z1, Real z2)
{
  return qag (z1, z2);
}

double RAD_OpticalDepthZ::integrand (double z)
{
  double u = uPZ (itsRADOD->getP (), z);
  double w = itsRADOD->getW (u);
  double mu = muPZ (itsRADOD->getP (), z);
  double wz = mu * w;
  double phi = itsRADOD->getPhi (wz);
  return u * u * phi / w;
}
