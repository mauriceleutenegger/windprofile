/***************************************************************************
    NumericalOpticalDepth.cpp   - Computes the optical depth numerically
                                by breaking the integral into pieces and
                                integrating the pieces along du or dz.

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
#include <iostream>

using namespace std;

const Real NumericalOpticalDepth::LARGE_OPTICAL_DEPTH = 1.e6;

const Real NumericalOpticalDepth::MU0 = 0.6;

const Real NumericalOpticalDepth::LARGE_P = 1.e4;

/*NumericalOpticalDepth::NumericalOpticalDepth 
(Real TauStar, Porosity* P, Velocity* V)
  : itsP (0.), itsTauStar (TauStar), isTransparent (false), 
    isPorous (false), itsPorosity (P), itsVelocity (V), 
    itsNumericalOpticalDepthZ (NULL), itsNumericalOpticalDepthU  (NULL)
    isHeII (false), nButterworth (3), uButterworth (5.)
{
  checkInput ();
  allocateNumericalOpticalDepthZU ();
  return;
  }*/

NumericalOpticalDepth::NumericalOpticalDepth 
(Real TauStar, bool HeII, Porosity* P, Velocity* V)
  : itsP (0.), itsTauStar (TauStar), isTransparent (false), 
    isPorous (false), itsPorosity (P), itsVelocity (V), 
    itsNumericalOpticalDepthZ (NULL), itsNumericalOpticalDepthU  (NULL),
    isHeII (HeII), nButterworth (3), uButterworth (0.2)
{
  /* enable this for debug:
  if (isHeII) {
    cout << "NumericalOpticalDepth: I am using HeII\n";
  } else {
    cout << "NumericalOpticalDepth: I am not using HeII\n";
  }
  */
  checkInput ();
  allocateNumericalOpticalDepthZU ();
  return;
}

NumericalOpticalDepth::~NumericalOpticalDepth ()
{
  freeNumericalOpticalDepthZU ();
  return;
}

void NumericalOpticalDepth::setParameters 
(Real TauStar, Porosity* P, Velocity* V)
{
  itsTauStar = TauStar;
  itsPorosity = P;
  itsVelocity = V;
  checkInput ();
}

void NumericalOpticalDepth::checkInput ()
{
  switch (compare (itsTauStar, 0.)) {
  case 1: isTransparent = false; break;
  case 0: isTransparent = true; break;
  default: cerr << "NumericalOpticalDepth: TauStar " << itsTauStar << "\n";
    isTransparent = true; break;
  }
  if (itsVelocity == 0) {
    cerr << "NumericalOpticalDepth: Velocity not initialized.\n";
    isTransparent = true;
    return;
  }
  if (compare (itsVelocity->getMinimumVelocity (), 0.) != 1) {
    cerr << "NumericalOpticalDepth: Velocity minimum not set.\n";
    isTransparent = true;
    return;
  }
  if (itsPorosity == 0) {
    isPorous = false;
  } else {
    isPorous = itsPorosity->getPorous ();
  }
  return;
}

void NumericalOpticalDepth::setTauStar (Real TauStar) 
{
  itsTauStar = TauStar;
  isTransparent = false;
  checkInput (); 
  return;
  // Note that the Porosity class will also be changed by the OpticalDepth class.
}

void NumericalOpticalDepth::allocateNumericalOpticalDepthZU ()
{
  itsNumericalOpticalDepthZ = new NumericalOpticalDepthZ (this);
  itsNumericalOpticalDepthU = new NumericalOpticalDepthU (this);
  return;
}

void NumericalOpticalDepth::freeNumericalOpticalDepthZU ()
{
  delete itsNumericalOpticalDepthZ;
  delete itsNumericalOpticalDepthU;
  return;
}

// This algorithm makes use of two integrals for the optical depth:
// One is in the standard dz coordinate, while the other is in the 
// inverse radial coordinate du. The du integral is useful for integrating 
// the outer region of the wind; a dz integral over this region does not 
// converge as rapidly. However, the du integral does poorly at small mu.
//
// Thus, the integral is broken up depending on mu of the starting coordinate.
// Furthermore, the process is more efficient if one makes use of the 
// symmetry of the integral at the midplane. Negative values of z are solved
// by recursively calling the function at |z| in combination with a call at
// z = 0.
Real NumericalOpticalDepth::getOpticalDepth (Real p, Real z)
{
  if (badCoordinates (p, z)) return LARGE_OPTICAL_DEPTH;
  if (isTransparent) return 0.;
  if (compare (z, 0.) == -1) {
    return (2. * getOpticalDepth (p, 0.) - getOpticalDepth (p, fabs(z)));
  }
  itsP = p;
  Real mu = muPZ (p, z);
  Real t = 0.;
  int status = 0;
  if (compare (mu, MU0) == 1) {
    t = itsNumericalOpticalDepthU->getOpticalDepth (p, z);
    status = itsNumericalOpticalDepthU->getStatus ();
    if (status) {
      cout << "NumericalOpticalDepth: NumericalOpticalDepthU returned status code " << status << "\n";
    }
  } else {
    if (compare (p, LARGE_P) == 1) {
      return itsTauStar * (M_PI_2 + atan (z / p)) / p;
    }
    t = itsNumericalOpticalDepthZ->getOpticalDepth (p,z);
    status = itsNumericalOpticalDepthZ->getStatus ();
    if (status) {
      cout << "NumericalOpticalDepth: NumericalOpticalDepthZ returned status code " << status << "\n";
    }
  }
  return itsTauStar * t;
}

// Approximate the He ionization fraction with a Butterworth filter
Real NumericalOpticalDepth::HeIIFilter (Real u) 
{
  Real answer = 1.;
  if (isHeII) { // Function for He ionization
    Real u2n = pow (u, 2*nButterworth); 
    Real ub2n = pow (uButterworth, 2*nButterworth);
    answer = 1. - sqrt (u2n / (u2n + ub2n)); 
  }
  return answer;
} 
