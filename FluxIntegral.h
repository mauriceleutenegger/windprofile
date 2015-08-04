/***************************************************************************
    FluxIntegral.h   - Integrates Lx dx
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

#ifndef MAL_FLUX_INTEGRAL_H
#define MAL_FLUX_INTEGRAL_H

#include <stdbool.h>
#include "xsTypes.h"
#include "mal_integration.h"
#include "Utilities.h"
#include "Lx.h"
#include "WindParameter.h"

class FluxIntegral : public Integral
{
 public:
  FluxIntegral (Lx* lx);
  ~FluxIntegral ();
  Real getFlux (Real x1, Real x2);
  Real getFlux (); // integrate over -1. < x < 1.
  double integrand (double x);
 private:
  Lx* itsLx;
  // To prevent copying or assignment;
  FluxIntegral (const FluxIntegral & I);
  FluxIntegral operator = (const FluxIntegral & I);
  Real itsXKink;
  Real itsXOcc;
};

#endif//MAL_FLUX_INTEGRAL_H
