/***************************************************************************
    NumericalOpticalDepthU.h   - Computes the optical depth numerically
                                 by doing the integral along du.

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

NumericalOpticalDepthU::NumericalOpticalDepthU (NumericalOpticalDepth* A)
  : Integral (), itsP (0.), itsNumericalOpticalDepth (A)
{
  return;
}

Real NumericalOpticalDepthU::getOpticalDepth (Real p, Real z)
{
  itsP = p;
  Real u = 1. / hypot (itsP, z);
  return qag (0., u);
}

double NumericalOpticalDepthU::integrand (double u)
{
  double f = itsNumericalOpticalDepth->HeIIFilter (u);
  double nu = itsP * u;
  double mu = sqrt (1. - nu * nu);
  double w = itsNumericalOpticalDepth->itsVelocity->getVelocity (u);
  if (!itsNumericalOpticalDepth->isPorous) return (f / (w * mu)); // smooth
  double PorosityFactor = 
    itsNumericalOpticalDepth->itsPorosity->getPorosityFactor (u, mu);
  return (f * PorosityFactor/ (w * mu));
}
