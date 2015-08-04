/***************************************************************************
    HeLikeRatio.cpp   - Calculates R = f/i as a function of radius.

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

#include "HeLikeRatio.h"
#include <iostream>

using namespace std;

HeLikeRatio::HeLikeRatio (Real R0, Real P, HeLikeType type)
  : itsR0 (R0), itsP (P), itsType (type)
{
  checkInput ();
  return;
}

void HeLikeRatio::setParameters (Real R0, Real P, HeLikeType type)
{
  itsR0 = R0; 
  itsP = P; 
  itsType = type;
  checkInput();
  return;
}

void HeLikeRatio::checkInput ()
{
  if (compare (itsR0, 0.) == -1) {
    cerr << "HeLikeRatio: Invalid R0 " << itsR0;
    itsR0 = 1.;
  }
  if (compare (itsP, 0.) == -1) {
    cerr << "HeLikeRatio: Invalid P " << itsP;
    itsP = 0.;
  }
  return;
}

Real HeLikeRatio::getHeLikeFactor (Real u)
{
  if (itsType == wResonance) return 1.;
  Real R = getR (u);
  if (itsType == yIntercombination) return (1. / (1. + R));
  if (itsType == zForbidden) return (R / (1. + R)); // forbidden
  cerr << "HeLikeRatio: Error in type " << itsType << "\n";
  return 0.;
  // x (important for high Z) is not currently supported - 060807 MAL
}

Real HeLikeRatio::getR0 ()
{
  return itsR0;
}

Real HeLikeRatio::getR (Real u)
{
  Real dilution = 1. - sqrt (1. - u * u); // actually dilution is twice this
  return (itsR0 / (1. + itsP * dilution)); // <-- but this is correct
}
