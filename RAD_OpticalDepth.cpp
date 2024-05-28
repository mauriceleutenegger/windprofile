/***************************************************************************
    RAD_OpticalDepth.cpp   - Calculates the optical depth to X-rays for a stellar
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

RAD_OpticalDepth::RAD_OpticalDepth (Velocity* V, Real DeltaE, Real Gamma,
				    Real Tau0, Real Vinfty) :
  Integral (),
  itsP (0.), itsZ0 (0.0), itsMu0 (0.0), itsU0 (0.0), itsW0 (0.0), itsWz0 (0.0),
  itsDeltaE (DeltaE), itsGamma (Gamma), itsTau0 (Tau0), itsVinfty (Vinfty),
  itsVelocity (V)

{
  // need to check how vinfty, gamma, deltaE are passed to this function
  return;
}

Real RAD_OpticalDepth::getOpticalDepth (Real p, Real z)
{
  initialize (p, z);
  return itsTau0 * qagiu (z);

}

void RAD_OpticalDepth::initialize (Real p, Real z)
{
  itsP = p;
  itsZ0 = z; // starting z for integral
  itsU0 = uPZ (itsP, itsZ0);
  itsMu0 = muPZ (itsP, itsZ0); // starting mu for integral
  //itsU0 = hypot (itsP, itsZ0);
  itsW0 = itsVelocity->getVelocity (itsU0);
  itsWz0 = itsW0 * itsMu0;
  return;
}

double RAD_OpticalDepth::integrand (double z)
{
  //double u = 1. / hypot (itsP, z);
  double u = uPZ (itsP, z);
  double w = itsVelocity->getVelocity (u);
  double mu = muPZ (itsP, z);
  double wz = mu * w;
  double deltaw = wz - itsWz0; //
  //double x = itsDeltaE + (1.0 - itsDeltaE) * (wz - itsWz0) * itsVinfty;
  double x = itsDeltaE + (1.0 - itsDeltaE) * deltaw * itsVinfty;
  double phi = gsl_ran_cauchy_pdf (x, itsGamma/2.);
  // width arg of cauchy is HWHM
  double answer = (u*u / w) * phi; // density * line profile
  // stuff - 1/density * phi
  return answer;
}

