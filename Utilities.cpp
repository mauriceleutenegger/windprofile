/***************************************************************************
    Utilities.cpp   - Utility functions for windprofile and Gaussians; 
      also, class Velocity.

                             -------------------
    begin				: December 2006
    copyright			: (C) 2006 by Maurice Leutenegger
    email				: maurice.a.leutenegger@nasa.gov
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

#include <iostream>
#include <float.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_const_cgsm.h>
#include "Utilities.h"

using namespace std;

static const Real CONST_HC_KEV_A = GSL_CONST_CGSM_PLANCKS_CONSTANT_H * 
    GSL_CONST_CGSM_SPEED_OF_LIGHT * 1.e5 / GSL_CONST_CGSM_ELECTRON_VOLT;

int compare (Real a, Real b)
{
  return gsl_fcmp (a, b, DBL_EPSILON);
}

Real muPZ (Real p, Real z)
{
  return (z / hypot (z, p));
}

bool isOcculted (Real p, Real z)
{
  Real r = hypot (p, z);
  bool impact =  (compare (p, 1.) != 1);
  bool behind = (compare (z, 0.) != 1);
  bool inside = (compare (r, 1.) != 1);
  return (impact && (behind || inside));
}

bool badCoordinates (Real p, Real z)
{
  if (compare (p, 0.) == -1){
    cerr << "Negative p = " << p << "not supported.\n";
    return true;
  }
  if (isOcculted (p, z)) {
    cerr << "Call with invalid coordinates is occulted;"
	 << "p = " << p << ", z = " << z << "\n";
    return true;
  }
  return false;
}

size_t BinarySearch (const RealArray& array, double value)
{
  size_t lower = 0;
  size_t upper = array.size ();
  while ((upper - lower) > 1) {
    size_t probe = (upper + lower) / 2;
    /*   if (compare (value, array[probe]) == 1) // value > array[probe] */
    if (value > array[probe])
      lower = probe;
    else
      upper = probe;
  }
  return lower;
}

bool isBackwards (const RealArray& energy)
{
  return (compare (energy[0], energy[1]) == 1);
}

bool isOutsideRange (const RealArray& energy, Real RestEnergy)
{
  bool isHigh (compare (RestEnergy, energy[energy.size () - 1]) == 1);
  bool isLow (compare (RestEnergy, energy[0]) == -1);
  return (isHigh || isLow);
}

/* Unit conversions: */


Real convert_KMS_C (Real x)
{
  return (x * 1.e5 / GSL_CONST_CGSM_SPEED_OF_LIGHT) ;
  // beta = v (km/s) * (1.e5 cm / km) / c (cm / s)
}

Real convert_A_keV (Real Wavelength_A)
{
  return CONST_HC_KEV_A / Wavelength_A;
}

Real ShiftWavelength 
(Real WavelengthA, Real DeltaVelocityKMS, Real DeltaWavelengthMA)
{
  Real DeltaVelocityC = convert_KMS_C (DeltaVelocityKMS);
  WavelengthA *= (1. + DeltaVelocityC);
  WavelengthA += DeltaWavelengthMA / 1000.;
  return WavelengthA;
}

/*--------------------------Velocity-------------------------------*/

Velocity::Velocity (Real beta, Real MinimumVelocity) 
  : itsBeta (beta), itsMinimumVelocity (MinimumVelocity)
{
  checkInput ();
  return;
}

void Velocity::setMinimumVelocity (Real MinimumVelocity)
{
  itsMinimumVelocity = MinimumVelocity;
  checkInput ();
  return;
}

void Velocity::checkInput ()
{
  if (compare (itsBeta, 0.) == -1) {
    cerr << "Velocity: invalid beta " << itsBeta << "; setting to 1.\n";
    itsBeta = 1.;
  }
  if (compare (itsMinimumVelocity, 0.) == -1) {
    cerr << "Velocity: invalid MinimumVelocity " << itsMinimumVelocity 
	 << "; setting to 0.\n";
    itsMinimumVelocity = 0.;
  }
  return;
}

Real Velocity::getVelocity (Real u)
{
  return (itsMinimumVelocity + (1. - itsMinimumVelocity) 
	  * pow (1. - u, itsBeta));
}
