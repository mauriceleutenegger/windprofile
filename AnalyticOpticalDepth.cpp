/***************************************************************************
    AnalyticOpticalDepth.cpp   - Computes wind optical depth analytically.

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

#include "AnalyticOpticalDepth.h"
#include <iostream>
#include <gsl/gsl_math.h>
#include "Utilities.h"

using namespace std;

AnalyticOpticalDepth::AnalyticOpticalDepth 
(Real TauStar, Real h, bool anisotropic, bool expansion)
  : itsTauStar (TauStar), itsTauClump0 (TauStar * h), isTransparent (false),
    isPorous (true), isAnisotropic (anisotropic), isExpansion (expansion)
{
  checkInput ();
}

void AnalyticOpticalDepth::checkInput ()
{
  switch (compare (itsTauStar, 0.)) {
  case 1: break;
  case 0: isTransparent = true; break;
  default: cerr << "OpticalDepth: Error in TauStar " << itsTauStar << "\n"; 
    isTransparent = true; break;
  }
  switch (compare (itsTauClump0, 0.)) {
  case 1: break;
  case 0: isPorous = false; break;
  default: cerr << "OpticalDepth: Error in TauClump0" << itsTauClump0 << "\n";
    isPorous = false; break;
  }
}

void AnalyticOpticalDepth::setParameters (Real TauStar, Real h)
{
  itsTauStar = TauStar;
  itsTauClump0 = TauStar * h;
  isTransparent = false;
  isPorous = true;
  checkInput ();
  return;
}

const Real AnalyticOpticalDepth::LARGE_OPTICAL_DEPTH = 1.e6;

Real AnalyticOpticalDepth::getOpticalDepth (Real p, Real z)
{
  if (isOcculted (p, z)) return LARGE_OPTICAL_DEPTH;
  if (isTransparent) return 0.;
  if (!isPorous) {
    return itsTauStar * Smooth (p, z);
  } else {
    if (isAnisotropic) {
      return itsTauStar * Anisotropic (p, z);
    } else { // Isotropic
      if (isExpansion) {
	return itsTauStar * Expansion (p, z);
      } else { // Stretch
	return itsTauStar * Isotropic (p, z);
      }
    }
  }
}

/**************************** Smooth Optical Depth ***********************/

Real AnalyticOpticalDepth::Smooth (Real p, Real z)
{
  const Real epsilon = 0.0001;
  mu = muPZ (p, z);
  zstar = sqrt (fabs (1. - p * p));
  //  Real argument = fabs (zstar / mu);
  //  if ((compare (argument, epsilon) == -1) && (compare (z, 0.) == 1)){
  Real argument = zstar / mu;
  if ((gsl_fcmp (argument, 0., epsilon) == 0) && (compare (z, 0.) == 1)) {
    SmoothA1 X (p, z);
    return X.sumSeries ();
  } else if (compare (p, 1.) == 1) {
    return SmoothG1 (p, z);
  } else {
    return SmoothL1 (p, z);
  }
}

inline Real AnalyticOpticalDepth::SmoothG1 (Real p, Real z) // p > 1
{
  return (M_PI_2 + atan (1. / zstar) - atan (z / zstar) 
	  - atan (mu / zstar)) / zstar;
}

Real AnalyticOpticalDepth::SmoothL1 (Real p, Real z) // 0 < p < 1
{
  Real r = hypot (p, z);
  Real a = r / (r * r - 1.);
  Real argument = (mu + zstar) * (z + zstar) * a / (1. + zstar);
  return (log (argument) / zstar);
}

/******************* Expansion Optical Depth Calculations ***************/

Real AnalyticOpticalDepth::Expansion (Real p, Real z)
{
  switch (compare (itsTauClump0, 1.)) {
  case 1:
    return ExpansionG1 (p, z);
    break;
  case 0:
    return Expansion1 (p, z);
    break;
  case -1:
    return ExpansionL1 (p, z);
  default:
    cerr << "OpticalDepth::Expansion: Error in TauClump0 = " << 
      itsTauClump0 << "\n";
    return 0.;
    break;
  }
}

/*--------------- Porous, expansion, TauClump0 > 1, very porous -----------*/
Real AnalyticOpticalDepth::ExpansionG1 (Real p, Real z)
{
  Real s = itsTauClump0 - 1.;
  Real mu = muPZ (p, z);
  Real zh;
  if (compare (p, 0.) == 0) return (log (1. + s / z) / s); // p = 0
  // Note: can s < 0? Is there a point to check if p = 0?
  switch (compare (p, s)) {
  case 1: // p > s
    zh = sqrt (p*p - s*s);
    return ((M_PI_2 + atan (-1. * s / zh) 
	  - atan (-1. * s * mu / zh) - atan (z / zh)) / zh);
    break;
  case 0: // p = s
    return (((mu - 1.) / (s * mu)) + (1. / z));
    break;
  case -1: // p < s
    zh = sqrt (s*s - p*p);
    if (compare (z, zh) == 0) { // z = zh, p < s
      return log (1 + zh / s) / zh;
    } else { // z != zh, p < s
      Real a = (s + zh) / (s - zh);
      Real b = (s * mu - zh) / (s * mu + zh);
      Real c = (z + zh) / (z - zh);
      return (log (a * b * c) / (2. * zh));
    }
    break;
  default:
    cerr << "OpticalDepth::ExpansionG1: Error in p = " << p << 
      "and s = " << s << "\n";
    return 0.;
    break;
  }
}

/*------------ Porous, expansion, TauClump0 = 1, marginally porous --------*/
inline Real AnalyticOpticalDepth::Expansion1 (Real p, Real z)
{
  if (compare (p, 0.) == 0) return (1. / z); // p = 0
  return  ((M_PI_2 - atan (z / p)) / p);
}

/*------------ Porous, expansion, TauClump0 < 1, not very porous ----------*/
Real AnalyticOpticalDepth::ExpansionL1 (Real p, Real z)
{
  Real s = 1. - itsTauClump0;
  Real zh, a, b, c;
  Real mu = muPZ (p, z);
  if (compare (p, 0.) == 0) return ((log(z / (z - s))) / s); // p = 0
  switch (compare (p, s)) {
  case 1: // p > s
    zh = sqrt (p*p - s*s);
    return ((M_PI_2 + atan (s / zh) - atan (s * mu / zh) - atan (z / zh)) / zh);
    break;
  case 0: // p = s
    return (((1. - mu) / (mu * s)) + (1. / z));
    break;
  case -1: // p < s
    zh = sqrt (s*s - p*p);
    a = (s - zh) / (s + zh);
    b = (s * mu + zh) / (s * mu - zh);
    c = (z + zh) / (z - zh);
    return (log (a * b * c) / (2. * zh));
    break;
  default:
    cerr << "OpticalDepth::ExpansionL1: Error in p = " << p << 
      "and s = " << s << "\n";
    return 0.;
    break;
  }
}

/****************** Porous, stretch, isotropic *********************/

Real AnalyticOpticalDepth::Isotropic (Real p, Real z)
{
  return ((Smooth (p, z) + Isotropic1 (p, z)) / (1. + itsTauClump0));
}

Real AnalyticOpticalDepth::Isotropic1 (Real p, Real z)
{
  const Real epsilon = 1.e-4;
  Real s2 = itsTauClump0;
  Real s = sqrt (s2);
  Real zh = hypot (p, s);
  Real r2 = p * p + z * z;
  mu = muPZ (p, z);
  Real a, y;
  if (gsl_fcmp (zh, 0., epsilon) != 1) {
    //  if (compare (zh, epsilon) != 1) {
    IsotropicSeries X (zh / z);
    a = (s2 / z) * X.sumSeries ();
  } else {
    a = (s2 / zh) * (M_PI_2 - atan (z / zh));
  }
  if (compare (p, s) == -1) {
    y = 1. / hypot (1., p / s);
  } else {
    y = s / zh;
  }
  Real b = y * atanh (zh * s / (s2 + r2 * (1. + mu)));
  return (a + b);
}

/****************** Porous, stretch, anisotropic ********************/

Real AnalyticOpticalDepth::Anisotropic (Real p, Real z)
{
  Real r1squared = (p * p + hypot (p * p, 2. * itsTauClump0)) / 2.;
  Real r1 = sqrt (r1squared);
  Real z1 = sqrt (r1squared - p * p);
  Real r = hypot (p, z);
  if (compare (z, z1) == 1) // z > z1
    return Smooth (p, z);
  if (compare (z, 0.) >= 0) // 0 <= z < z1
    return (Smooth (p, z1) + Anisotropic1 (r, r1) / itsTauClump0);
  if (compare (z, -1. * z1) >= 0) // -z1 <= z < 0
    return (Smooth (p, z1) + (2. * Anisotropic1 (p, r1) 
			      - Anisotropic1 (r, r1)) / itsTauClump0);
  else // z < -z1
    return (2. * (Smooth (p, z1) + Anisotropic1 (p, r1) / itsTauClump0) 
	    - Smooth (p, fabs(z)));
}

inline Real AnalyticOpticalDepth::Anisotropic1 (Real ra, Real rb)
{
  return (rb - ra + log ((rb - 1.) / (ra - 1.)));
}

