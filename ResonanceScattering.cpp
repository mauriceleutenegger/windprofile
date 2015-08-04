/***************************************************************************
    ResonanceScattering.cpp   - Calculates the effects of resonance scattering 
                              on a stellar wind X-ray emission line profile.

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

#include "ResonanceScattering.h"
#include <iostream>
#include <complex>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_exp.h>

using namespace std;

ResonanceScattering::ResonanceScattering
(Real Tau0Star, Real BetaSobolev, bool OpticallyThick, Velocity* V) 
  : itsTau0Star (Tau0Star), itsBetaSobolev (BetaSobolev),
    isOpticallyThick (OpticallyThick), isOpticallyThin (true), itsVelocity (V)
{
  checkInput ();
  return;
}

void ResonanceScattering::setParameters 
(Real Tau0Star, Real BetaSobolev, bool OpticallyThick, Velocity* V) 
{
  itsTau0Star = Tau0Star; 
  itsBetaSobolev = BetaSobolev; 
  isOpticallyThick = OpticallyThick; 
  itsVelocity = V; 
  checkInput ();
  return;
}

void ResonanceScattering::checkInput ()
{
  isOpticallyThin = true;
  switch (compare (itsTau0Star, 0.)) {
  case 1: isOpticallyThin = false; break;
  case 0: break;
  default: cerr << "ResonanceScattering: Bad Tau0Star " << itsTau0Star << "\n";
    break;
  }
  if (compare (itsBetaSobolev, 0.) == -1) {
    cerr << "ResonanceScattering: Bad BetaSobolev " << itsBetaSobolev << "\n";
    itsBetaSobolev = 0.;
  }
  if (isOpticallyThick) isOpticallyThin = false;
  return;
}

Real ResonanceScattering::getEscapeProbability (Real u, Real mu)
{
  if (isOpticallyThin) return 1.;
  Real sigma = (itsBetaSobolev * u) / (1. - u) - 1.;
  Real AngleFactor = (1. + sigma * mu * mu);
  if (isOpticallyThick) {
    return (AngleFactor / (1. + sigma / 3.));
  }
  Real w = itsVelocity->getVelocity (u);
  Real Tau0 = itsTau0Star * u / (w * w);
  if (compare (AngleFactor, 0.) != 1) {
    return 0.;
  }
  Real TauMu = -1. * Tau0 / AngleFactor;
  Real p = gsl_sf_exprel (TauMu);
  //  Real p = expm1 (TauMu) / TauMu; // 1 - e^-t / t
  Real PAverage = getPAverage (Tau0, sigma);
  return (p / PAverage);
}

/* Castor, Radiation hydrodynamics, pp 128-129. */
Real ResonanceScattering::getPAverage (Real Tau0, Real sigma)
{
  complex<Real> z (-0.97515, -1.193464);
  complex<Real> rbz (-0.5, 0.01193391);
  // if sigma = 0, analytically simplify first:
  if (compare (sigma, 0.) == 0) {
    complex<Real> temp = rbz / ((Tau0 / z) - 1.);
    return ((Real) (2. * real (temp)));
  }
  // could do sigma = 0 and directly calculate integral,
  // but that might be bad for continuity of result with non-zero sigma
  complex<Real> t = sqrt (((Tau0 / z) - 1.) / sigma);
  complex<Real> argument = (t - 1.) / (t + 1.);
  complex<Real> tlog = Tau0 * log (argument);
  complex<Real> temp = rbz * ((tlog / (t * z * (sigma * 2.))) + 1.);
  return ((Real) (-2. * real (temp)));
}
