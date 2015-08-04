/***************************************************************************
    IntegratedLuminosity.h
       integrates total luminosity over inverse radius 
       from 0 to u0
     

                             -------------------
    begin				: July 2008
    copyright			: (C) 2008 by Maurice Leutenegger
    email				: Maurice.A.Leutenegger@nasa.gov
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

#ifndef MAL_INTEGRATED_LUMINOSITY_H
#define MAL_INTEGRATED_LUMINOSITY_H

#include <gsl/gsl_math.h>
#include "xsTypes.h"
#include "Utilities.h"
#include "AngleAveragedTransmission.h"
//#include "OpticalDepth.h"
//#include "ResonanceScattering.h"
//#include "HeLikeRatio.h"
#include "mal_integration.h"
//#include "UxRoot.h"


class IntegratedLuminosity : public Integral {
 public:
  IntegratedLuminosity 
    (Real q, Real itsU0, Real itsUmin, Velocity* V, 
     AngleAveragedTransmission* T);
  ~IntegratedLuminosity ();
  Real getLuminosity ();
  double integrand (double u);
  void setTransparentCore 
    (bool TransparentCore) {isTransparentCore=TransparentCore;}
  void setU0 (double u0) {itsU0 = u0;}
 private:
  bool isTransparentCore;
  Real itsQ;
  Real itsU0;
  Real itsUmin;
  Velocity* itsVelocity;
  AngleAveragedTransmission* itsAngleAveragedTransmission;
};

#endif//MAL_INTEGRATED_LUMINOSITY_H
