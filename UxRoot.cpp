/***************************************************************************
    UxRoot.cpp   - Finds the u coordinate at which p = 1 + epsilon for a given
                 value of x. This allows the avoidance of the occulted part
                 of the Lx integral.

                             -------------------
    begin				: March 2007
    copyright			: (C) 2007 by Maurice Leutenegger
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

#include "UxRoot.h"
#include <iostream>

using namespace std;

static const Real epsilon = 1.e-5;
static const Real p = 1. + epsilon;
static const Real p2 = p * p;

UxRoot::UxRoot (Velocity* V, Real  x) :
  RootFinderNewton (0., 1.e-3, 0.), itsVelocity (V), itsX2 (0)
{
  setX (x);
  return;
}

double UxRoot::f (double u)
{
  Real w = itsVelocity->getVelocity (u);
  return w * w * (1. - p2 * u * u) - itsX2;
}

double UxRoot::df (double u)
{
  Real w = itsVelocity->getVelocity (u);
  Real beta = itsVelocity->getBeta ();
  return -2. * w * w * (beta + p2 * u - (1. + beta) * u * u * p2) / (1. - u);
}

void UxRoot::fdf (double u, double& y, double& dy)
{
  Real w = itsVelocity->getVelocity (u);
  Real beta = itsVelocity->getBeta ();
  Real w2 = w * w;
  Real u2 = u * u;
  y = w2 * (1. - p2 * u2) - itsX2;
  dy = -2. * w2 * (beta + p2 * u - (1. + beta) * u2 * p2) / (1. - u);
  return;
}

/*
double UxRoot::f (double u)
{
  Real w = itsVelocity->getVelocity (u);
  return w * w * (1. - u * u) - itsX2;
}

double UxRoot::df (double u)
{
  Real w = itsVelocity->getVelocity (u);
  Real beta = itsVelocity->getBeta ();
  return -2. * w * w * (beta + (1. + beta) * u);
}

void UxRoot::fdf (double u, double& y, double& dy)
{
  Real w = itsVelocity->getVelocity (u);
  Real beta = itsVelocity->getBeta ();
  y = w * w * (1. - u * u) - itsX2;
  dy = -2. * w * w * (beta + (1. + beta) * u);
  return;
}
*/

void UxRoot::setX (Real x)
{
  itsX2 = x * x;
  Real beta = itsVelocity->getBeta ();
  setGuess (1. - pow (fabs (x), 1. / beta));
  return;
}
