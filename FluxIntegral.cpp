/***************************************************************************
    FluxIntegral.cpp   - Integrates Lx dx
                       L (x1, x2) = int^x2_x1 dx Lx

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

#include "FluxIntegral.h"
#include <iostream>

using namespace std;

FluxIntegral::FluxIntegral (Lx* lx)
  : Integral (), itsLx (lx), itsXKink (itsLx->getXKink ()), 
    itsXOcc (itsLx->getXOcc ())
{
  return;
}

FluxIntegral::~FluxIntegral ()
{
  return;
}

double FluxIntegral::integrand (double x)
{
  return itsLx->getLx (x);
}

/* 
   Integrates Lx over [y,x].
   Will mostly be called with y < x, but if not, it swaps them.
   Note: is there code somewhere else to deal with an unresolved profile?
   (The code below only deals with the case where it's known to be 
   completely unresolved in advance.)
*/
Real FluxIntegral::getFlux (Real x, Real y)
{
  if (compare (y, x) == 1) {
    Real temp = x;
    x = y;
    y = temp;
  }
  bool xInRange = (compare (fabs(x), 1.) == -1);
  bool yInRange = (compare (fabs(y), 1.) == -1);
  if (!xInRange && !yInRange) return 0.;
  if (!xInRange) x = 1.;
  if (!yInRange) y = -1.;
  if ((compare (itsXKink, y) == 1) && (compare (itsXKink, x) == -1)) {
    return (qag (y, itsXKink) + qag (itsXKink, x));
  }
  if ((compare (itsXOcc, y) == 1) && (compare (itsXOcc, x) == -1)) {
    return (qag (y, itsXOcc) + qag (itsXOcc, x));
  }
  return qag (y, x);
  /*
    The kink coordinate gives the point at which many profiles have a kink at 
    negative x on the blue side of the profile. This kink occurs at the x 
    where the cutoff u0 starts to become important. (At zero optical depth, 
    the profile becomes flat at this point).
    The "occ" coordinate gives the point in x where occultation begins to 
    become important. I choose to place this not at p = 1 but at p = 1 + e, 
    to avoid edge effects.
   */
  /*
  if (xInRange && yInRange) {
    return qag (y, x);
  } else if (xInRange && !yInRange) {
    return qag (-1., x);
  } else if (!xInRange && yInRange) {
    return qag (y, 1.);
  } else {
    return 0.;
  }
  */
}

// Integrate on x in [-1:1]
Real FluxIntegral::getFlux ()
{
  return (qag (-1., itsXKink) + qag (itsXKink, itsXOcc) + qag (itsXOcc, 1.));
}
