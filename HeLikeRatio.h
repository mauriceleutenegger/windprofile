/***************************************************************************
    HeLikeRatio.h   - Calculates R = f/i as a function of radius.

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

#ifndef HE_LIKE_RATIO_H
#define HE_LIKE_RATIO_H

#include <gsl/gsl_math.h>
#include "xsTypes.h"
#include "Utilities.h"

class HeLikeRatio
{
 public:
  HeLikeRatio (Real R0 = 1., Real P = 0., HeLikeType type = wResonance);
  void setParameters (Real R0, Real P, HeLikeType type);
  void setHeLikeType (HeLikeType type) {itsType = type; return;}
  Real getHeLikeFactor (Real u);
  Real getR0 ();
 private:
  Real itsR0; // ratio with no photoexcitation
  Real itsP; // phi_* / phi_c
  HeLikeType itsType; // type of line
  Real getR (Real u); // ratio as a function of radius
  void checkInput ();
};

#endif//HE_LIKE_RATIO_H
