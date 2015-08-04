/***************************************************************************
    Porosity.cpp   - Computes the porosity factor P = kappa_eff / kappa.

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

#include "Porosity.h"
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_exp.h>

using namespace std;

Porosity::Porosity 
(Real TauClump0, bool Anisotropic, bool Rosseland)
  : itsTauClump0 (TauClump0), isAnisotropic (Anisotropic), 
    isRosseland (Rosseland), isPorous (false)
{
  checkInput ();
}

void Porosity::setParameters 
(Real TauClump0, bool Anisotropic, bool Rosseland) 
{
  itsTauClump0 = TauClump0;
  isAnisotropic = Anisotropic;
  isRosseland = Rosseland; 
  isPorous = false;
  checkInput ();
  return;
}

void Porosity::setParameters (Real TauClump0) 
{
  itsTauClump0 = TauClump0;
  isPorous = false;
  checkInput ();
  return;
}

void Porosity::checkInput ()
{
  switch (compare (itsTauClump0, 0.)) {
  case 1:
    isPorous = true; return; break;    
  case 0:
    isPorous = false; return; break;
  default:
    cerr << "Porosity: TauClump = " << itsTauClump0 << "; setting to 0.\n";
    itsTauClump0 = 0.; isPorous = false; isAnisotropic = false; 
    isRosseland = false; return; break;
  }
  return;
}

Real Porosity::getTauClump (Real u) 
{
  return itsTauClump0 * u * u; // isotropic
} 

Real Porosity::getPorosityFactor (Real u, Real mu)
{
  mu = fabs (mu);
  if (!isPorous) return 1.;
  Real TauClump = getTauClump (u); // isotropic clump optical depth
  if (isRosseland) {
    if (isAnisotropic) {
      return (mu / (mu + TauClump));
    } else { //istropic
      return (1. / (1. + TauClump));
    }
  } else { //regular bridging
    if (isAnisotropic) {
      // in case mu is small:
      if (compare (mu / TauClump, 1.e-10) != 1) return (mu / TauClump);
      TauClump /= mu;
    } 
    return gsl_sf_exprel (-1. * TauClump); // (1 - e^-t) / t
  }
}
