/***************************************************************************
    Lx.h   - Calculates Lx for a stellar wind X-ray emission line profile.

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

#ifndef MAL_LX_H
#define MAL_LX_H

#include <gsl/gsl_math.h>
#include "xsTypes.h"
#include "Utilities.h"
#include "OpticalDepth.h"
#include "ResonanceScattering.h"
#include "HeLikeRatio.h"
#include "mal_integration.h"
#include "UxRoot.h"

class Lx : public Integral
{
 public:
  Lx (Real q, Real U0, Real Umin, Real beta, HeLikeType type, Velocity* V, 
      HeLikeRatio* He, ResonanceScattering* RS, OpticalDepth* Tau);
  Lx (Real q, Real U0, Real Umin, Real beta, Real kappaRatio, HeLikeType type, 
      Velocity* V, HeLikeRatio* He, ResonanceScattering* RS, 
      OpticalDepth* Tau, OpticalDepth* TauHeII);
  ~Lx ();
  void setHeLikeType (HeLikeType type);
  void setTransparent () {isTransparent = true; return;}
  void notTransparent () {isTransparent = false; return;}
  Real getLx (Real x);
  Real getXKink ();
  Real getXOcc ();
  double integrand (double u);
 private:
  Real itsX;
  Real itsQ;
  Real itsU0;
  Real itsUmin;
  Real itsBeta;
  Real itsKappaRatio;
  bool isTransparent;
  bool isHeLike;
  bool isHeII;
  HeLikeType itsHeLikeType;
  Velocity* itsVelocity;
  HeLikeRatio* itsHeLikeRatio;
  ResonanceScattering* itsResonanceScattering;
  OpticalDepth* itsOpticalDepth;
  OpticalDepth* itsOpticalDepthHeII;
  UxRoot* itsUxRoot;
  double integrand0 ();
  void checkInput ();
  void allocateClasses ();
  void freeClasses ();
  // To prevent copying and assignment:
  Lx (const Lx & l);
  Lx operator = (const Lx & l);
};

#endif//MAL_LX_H
