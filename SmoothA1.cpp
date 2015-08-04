/***************************************************************************
    SmoothA1.h   - Calculates part of the analytic expression for
                   smooth optical depth using a series expansion.

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

#include "SmoothA1.h"
#include "Utilities.h"
#include <iostream>

using namespace std;

SmoothA1::SmoothA1 (double pp, double zz) 
  : Series (), p (pp), z (zz), mu2 (1.), muTerm (1.), z2 (z * z),
    zTerm (1.), zStar2 (1. - p * p), zStarTerm (1.)
{
  checkInput ();
  initialize ();
  return;
}

void SmoothA1::checkInput () 
{
  const double epsilon = 0.1;
  if (compare (z, 0.) != 1) {
    cout << "SmoothA1: bad z " << z << "\n";
    z = 1.;
  }
  if (compare (fabs (p - 1.), epsilon) != -1) {
    cout << "SmoothA1: bad p " << p << "\n";
    p = 1.;
  }
  return;
}

void SmoothA1::initialize ()
{
  double mu = muPZ (p, z);
  mu2 = mu * mu;
  muTerm = 1. / mu;
  zTerm = 1. / z;
  return;
}

double SmoothA1::getTerm () 
{
  return (zStarTerm / (2. * getN () + 1.)) * (zTerm + muTerm - 1.);
}
// need to check this - doesn't n in class Series start at one?
// but it needs to start at zero!

void SmoothA1::iterate ()
{
  zStarTerm *= zStar2;
  muTerm /= mu2;
  zTerm /= z2;
  return;
}
