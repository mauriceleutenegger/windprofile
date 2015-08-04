/***************************************************************************
    Utilities.h   - Utility functions for windprofile and Gaussians; 
      also, class Velocity.

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

#ifndef MAL_WINDPROF_UTILITIES_H
#define MAL_WINDPROF_UTILITIES_H

#include <stdbool.h>
#include "xsTypes.h"

enum HLikeType {single, blue, red};
//enum HeLikeType {resonance, intercombination, forbidden};
enum HeLikeType {wResonance, xIntercombination, yIntercombination, zForbidden};

int compare (Real a, Real b);
Real muPZ (Real p, Real z);
bool isOcculted (Real p, Real z);
bool badCoordinates (Real p, Real z);
size_t BinarySearch (const RealArray& array, Real value);
bool isBackwards (const RealArray& energy);
bool isOutsideRange (const RealArray& energy, Real RestEnergy);

/* Unit conversions: */

Real convert_KMS_C (Real x);
Real convert_A_keV (Real Wavelength_A);

Real ShiftWavelength 
(Real WavelengthA, Real DeltaVelocityKMS, Real DeltaWavelengthMA);


class Velocity
{
 public:
  Velocity (Real beta = 1., Real MinimumVelocity = 0.);
  void setBeta (Real beta) {itsBeta = beta; checkInput (); return;}
  void setMinimumVelocity (Real MinimumVelocity);
  Real getVelocity (Real u);
  Real getMinimumVelocity () const {return itsMinimumVelocity;}
  Real getBeta () const {return itsBeta;}
 private:
  Real itsBeta;
  Real itsMinimumVelocity;
  void checkInput ();
};

#endif//MAL_WINDPROF_UTILITIES_H
