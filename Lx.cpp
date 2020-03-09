/***************************************************************************
    Lx.cpp   - Calculates Lx for a stellar wind X-ray emission line profile.

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

#include "Lx.h"
#include <iostream>

using namespace std;

Lx::Lx (Real q, Real U0, Real Umin, Real beta, HeLikeType type, 
	Velocity* V, HeLikeRatio* He, ResonanceScattering* RS, OpticalDepth* Tau)
  : Integral (), itsQ (q), itsU0 (U0), itsUmin (Umin), itsBeta (beta), itsKappaRatio (0.),
    isTransparent (false), isHeLike (false), isHeII (false),
    itsHeLikeType (type), itsVelocity (V), itsHeLikeRatio (He), 
    itsResonanceScattering (RS), itsOpticalDepth (Tau), 
    itsOpticalDepthHeII (NULL), itsUxRoot (0) 
{
  checkInput ();
  allocateClasses ();
  return;
}

Lx::Lx (Real q, Real U0, Real Umin, Real beta, Real kappaRatio, HeLikeType type, 
	Velocity* V, HeLikeRatio* He, ResonanceScattering* RS, 
	OpticalDepth* Tau, OpticalDepth* TauHeII)
  : Integral (), itsQ (q), itsU0 (U0), itsUmin (Umin), itsBeta (beta), itsKappaRatio (kappaRatio),
    isTransparent (false), isHeLike (false), isHeII (true),
    itsHeLikeType (type), itsVelocity (V), itsHeLikeRatio (He), 
    itsResonanceScattering (RS), itsOpticalDepth (Tau), 
    itsOpticalDepthHeII (TauHeII), itsUxRoot (0) 
{
  checkInput ();
  allocateClasses ();
  return;
}

void Lx::checkInput ()
{
  if (compare (itsQ, -1.) != 1) {
    cerr << "Lx: invalid q " << itsQ << "\n";
    itsQ = 0.;
  }
  if (compare (itsU0, 0.) != 1) {
    cerr << "Lx: invalid U0 " << itsU0 << "\n";
    itsU0 = 0.5;
  }
  if (compare (itsBeta, 0.) == -1) {
    cerr << "Lx: invalid beta " << itsBeta << "\n";
    itsBeta = 1.;
  }
  if (itsHeLikeRatio != 0) isHeLike = true;
  return;
}

void Lx::allocateClasses ()
{
  itsUxRoot = new UxRoot (itsVelocity, 0.5);
}

Lx::~Lx ()
{
  freeClasses ();
  return;
}

void Lx::freeClasses ()
{
  delete itsUxRoot;
}

// Unfortunately, Lx needs to know this so that it can decide whether
// to use resonance scattering or He-like ratios. 
// Is there a way to remove the necessity for it to know this?
void Lx::setHeLikeType (HeLikeType type)
{
  itsHeLikeType = type;
  return;
}

// Is there a different MIN function available? 
// Perhaps in one of the c++ libraries.
Real Lx::getLx (Real x)
{
  if (compare (fabs (x), 1.) != -1) return 0.;
  itsX = x;
  double Ux = 1. - pow (fabs(itsX), 1. / itsBeta);
  Ux = GSL_MIN_DBL (Ux, itsU0);
  if (compare (x, getXOcc ()) != -1) {
    itsUxRoot->setX (x);
    double UxOcc = itsUxRoot->findRoot ();
    Ux = GSL_MIN_DBL (Ux, UxOcc);
    // This is to make sure we start the integral at p = 1 in the 
    // red hemisphere.
  }
  /* If Umin > 0, need to check that Ux > Umin;
     otherwise there is no emission. */
  if (compare (Ux, itsUmin) != 1) {
    return 0.;
  }
  double answer = 0.;
  /* If Umin > 0., ignore convergence issues (see below), 
     since the Umin cutoff will likely solve the problem. */
  if (compare (itsUmin, 0.) == 1) {
    answer = qagp (itsUmin, Ux);
    return answer;
  }
  /* Split the integral for q < -0.5, since it seems to have trouble 
     converging sometimes (esp. for small tau_* and small x).
     The convergence problem occurs as u->0, so a small part of that end
     of the integral is split off. */
  if (compare (itsQ, -0.5) == -1) {
    double answer1 = qagp (Ux/10., Ux);
    double answer2 = qagp (0., Ux/10.);
    answer = answer1 + answer2;
  } else {
    answer = qagp (0., Ux);
  }
  return answer;
  //return qagp (0., Ux);
  // it's better to avoid the ambiguity of what happens at u = 0.
  /*
  // If q < 0, use qagp to handle u = 0 singularity.
  if (compare (itsQ, 0.) == -1) {
    return qagp (0., Ux);
  } else {
    return qag (0., Ux);
  }
  */
}

Real Lx::getXKink ()
{
  return -1. * pow (1. - itsU0, itsBeta);
}

Real Lx::getXOcc ()
{
  return itsVelocity->getVelocity (itsU0) * sqrt (1. - itsU0 * itsU0);
}

double Lx::integrand (double u)
{
  if (compare (u, 0.) == 0) return integrand0 ();
  Real w = itsVelocity->getVelocity (u);
  Real mu = -1. * itsX / w; 
  Real p = sqrt (1. - mu * mu) / u;
  Real z = mu / u;
  if  (isOcculted (p, z)) return 0.;
  // continuum optical depth
  Real Transmission = 1.;
  if (!isTransparent) { // Transparent is a flag to allow for easy
    // calculation of a profile with zero optical depth
    Real tau = itsOpticalDepth->getOpticalDepth (p, z);
    if (isHeII) { // opacity from He++ recombining to He+
      tau += itsKappaRatio * itsOpticalDepthHeII->getOpticalDepth (p, z);
    }
    Transmission = exp (-1. * tau);
  }
  /*
    This should be recoded so that HeLikeType is checked only for a He-like
    line (even though it works fine to just set type=wResonance).
  */
  Real HeLikeFactor = 1.;
  if (isHeLike && (itsHeLikeType != wResonance)) {
    HeLikeFactor = itsHeLikeRatio->getHeLikeFactor (u, itsHeLikeType);
  }
  Real EscapeProbability = 1.;
  /*
    This should also be recoded so that it is always eligible for RS unless
    it's specifically a non-resonance He-like line.
   */
  if (!isHeLike || (itsHeLikeType == wResonance)) { 
    EscapeProbability = itsResonanceScattering->getEscapeProbability (u, mu);
  }
  double Integrand = pow (u, itsQ) / gsl_pow_3 (w); // emission
  Integrand *= Transmission; // absorption
  Integrand *= HeLikeFactor; // radial dependence of f/i ratio
  Integrand *= EscapeProbability; // resonance scattering
  return Integrand;
}

// Integrand for u = 0.
double Lx::integrand0 ()
{
  switch (compare (itsQ, 0.)) {
  case 1: // q > 0.
    return 0.; // 0.^q = 0.
    break;
  case 0: // q = 0.
    return 1.; // tau = 0., w = 1., 0.^0. = 1.
    break;
  case -1: // q < 0, NAN; shouldn't be called.
  default:
    cerr << "LxIntegral::Integrand0: Error in q = " << itsQ << "\n";
    return 1.e10;
    break;
  }
}



