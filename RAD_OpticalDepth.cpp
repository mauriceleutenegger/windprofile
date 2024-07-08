/***************************************************************************
    RAD_OpticalDepth.h   - Calculates the optical depth to X-rays from 
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

const Real RAD_OpticalDepth::LARGE_OPTICAL_DEPTH = 1.e6;

const Real RAD_OpticalDepth::MU0 = 0.6;

RAD_OpticalDepth::RAD_OpticalDepth (Velocity* V, Real DeltaE, Real Gamma,
				    Real Tau0, Real Vinfty) :
  itsP (0.), itsZ0 (0.0), itsMu0 (0.0), itsU0 (0.0), itsW0 (0.0), itsWz0 (0.0),
  itsDeltaE (DeltaE), itsGamma (Gamma), itsTau0 (Tau0), itsVinfty (Vinfty),
  isTransparent (false),
  itsVelocity (V), itsRADODZ (NULL), itsRADODU (NULL)

{
  // need to check how vinfty, gamma, deltaE are passed to this function
  checkInput ();
  allocateRAD_OpticalDepthZU ();
  return;
}

void RAD_OpticalDepth::checkInput ()
{
  switch (compare (itsTau0, 0.)) {
  case 1: isTransparent = false; break;
  case 0: isTransparent = true; break;
  default: cerr << "RAD_OpticalDepth: Tau0 " << itsTau0 << "\n";
    isTransparent = true; break;
  }
  if (itsVelocity == 0) {
    cerr << "RAD_OpticalDepth: Velocity not initialized.\n";
    isTransparent = true;
    return;
  }
}

void RAD_OpticalDepth::allocateRAD_OpticalDepthZU ()
{
  itsRADODZ = new RAD_OpticalDepthZ (this);
  itsRADODU = new RAD_OpticalDepthU (this);
  return;
}

void RAD_OpticalDepth::freeRAD_OpticalDepthZU ()
{
  delete itsRADODZ;
  delete itsRADODU;
  return;
}

Real RAD_OpticalDepth::getOpticalDepth (Real p, Real z)
{
  if (badCoordinates (p, z)) return LARGE_OPTICAL_DEPTH;
  if (isTransparent) return 0.;

  initialize (p, z);

  Real t = 0.;

  // use U integral for mu0 > MU0
  if (compare (itsMu0, MU0) == 1) {
    t = itsRADODU->getOpticalDepth (z);
    return itsTau0 * t;
  }
  // z coordinates at MU0 given P
  Real zprime = MU0 * fabs (itsP) / sqrt (1. - MU0 * MU0);
  /* if mu0 < -MU0,
     break into U integral up to -MU0;
     Z integral from -MU0 to MU0;
     U integral from MU0 to 1 */
  if (compare (itsMu0, -1. * MU0) == -1) {
    t = itsRADODU->getOpticalDepth (z, -1. * zprime);
    t += itsRADODZ->getOpticalDepth (-1. * zprime, zprime);
    t += itsRADODU->getOpticalDepth (zprime);
    return itsTau0 * t;
  }
  // otherwise Z integral to MU0, then U integral
  t = itsRADODZ->getOpticalDepth (z, zprime);
  t += itsRADODU->getOpticalDepth (zprime);
  return itsTau0 * t;
}

void RAD_OpticalDepth::initialize (Real p, Real z)
{
  itsP = p;
  itsZ0 = z; // starting z for integral
  itsU0 = uPZ (itsP, itsZ0);
  itsMu0 = muPZ (itsP, itsZ0); // starting mu for integral
  itsW0 = itsVelocity->getVelocity (itsU0);
  itsWz0 = itsW0 * itsMu0;
  return;
}

// Lorentzian line profile used in both integrals
Real RAD_OpticalDepth::getPhi (Real wz)
{
  Real deltaw = wz - itsWz0;
  Real x = itsDeltaE + (1.0 - itsDeltaE) * deltaw * itsVinfty;
  Real phi = gsl_ran_cauchy_pdf (x, itsGamma / 2.);
  return phi;
}

/* allows exploring the behavior of the integrand in the dz integration
   using the python interface */
double RAD_OpticalDepth::getZIntegrand (double z)
{
  return itsRADODZ->integrand (z);
}
