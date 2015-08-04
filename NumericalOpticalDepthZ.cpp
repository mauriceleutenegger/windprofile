/***************************************************************************
    NumericalOpticalDepthZ.h   - Computes the optical depth numerically
                                 by doing the integral along dz.

                             -------------------
    begin				: December 2006
    copyright			: (C) 2006 by Maurice Leutenegger
    email				: maurice@astro.columbia.edu
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

#include "NumericalOpticalDepth.h"

NumericalOpticalDepthZ::NumericalOpticalDepthZ (NumericalOpticalDepth* A)
  : Integral (), itsP (0.), itsNumericalOpticalDepth (A)
{
  return;
}

Real NumericalOpticalDepthZ::getOpticalDepth (Real p, Real z1, Real z2)
{
  itsP = p;
  return qag (z1, z2);
}

Real NumericalOpticalDepthZ::getOpticalDepth (Real p, Real z)
{
  itsP = p;
  return qagiu (z);
}

Real NumericalOpticalDepthZ::getOpticalDepth (Real p)
{
  itsP = p;
  return qagi ();
}

double NumericalOpticalDepthZ::integrand (double z)
{
  double u = 1. / hypot (itsP, z);
  double f = itsNumericalOpticalDepth->HeIIFilter (u);
  double mu = z * u;
  double w = itsNumericalOpticalDepth->itsVelocity->getVelocity (u);
  if (!itsNumericalOpticalDepth->isPorous) return (f * u * u / w); // smooth
  double PorosityFactor = 
    itsNumericalOpticalDepth->itsPorosity->getPorosityFactor (u, mu);
  return (f * u * u * PorosityFactor / w);
}
