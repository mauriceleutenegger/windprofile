/***************************************************************************
    RAD_OpticalDepth.h   - Calculates the optical depth to X-rays for a stellar
                       wind along a ray p from point z.

                             -------------------
    begin				: May 2024
    copyright			: (C) 2024 by Gabe Grell and Maurice Leutenegger
    email				: gabriel.j.grell@nasa.gov maurice.a.leutenegger@nasa.gov
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

#ifndef RAD_OPTICAL_DEPTH_H
#define RAD_OPTICAL_DEPTH_H

#include <stdbool.h>
#include "xsTypes.h"
#include "Utilities.h"
#include "mal_integration.h"

class RAD_OpticalDepth : public Integral {
 public:
  RAD_OpticalDepth (Velocity* V, Real DeltaE, Real Gamma, Real Tau0,
		    Real Vinfty); // constructor
  //~RAD_OpticalDepth (); // destructor
  Real getOpticalDepth (Real p, Real z);
  void initialize (Real p, Real z); // setup to manually call integrand
  // initialize is also called by getOpticalDepth
  double integrand (double z);
 private:
  // add the characteristic optical depth as a parameter
  // use it in getOpticalDepth
  // get it from the constructor
  Real itsP;
  Real itsZ0;
  Real itsMu0;
  Real itsU0;
  Real itsW0;
  Real itsWz0;
  Real itsDeltaE; // in units of rest energy of absorber; (Eabs - Eem) / Eabs
  Real itsGamma; // same units as DeltaE; i.e. gamma / Eabs
  Real itsTau0;
  Real itsVinfty; // in units of speed of light
  Velocity* itsVelocity;
  // are allocators and deallocators needed?

  // To prevent copying and assignment:
  //RAD_OpticalDepth (RAD_OpticalDepth const & T);
  //RAD_OpticalDepth operator = (RAD_OpticalDepth const & T);

};

# endif//RAD_OPTICAL_DEPTH_H
