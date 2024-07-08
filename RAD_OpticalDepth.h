/***************************************************************************
    RAD_OpticalDepth.h   - Calculates the optical depth to X-rays from 
                           Resonant Auger Destruction (RAD) for a stellar
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

class RAD_OpticalDepth;

class RAD_OpticalDepthZ : public Integral {
public:
  RAD_OpticalDepthZ (RAD_OpticalDepth* A);
  ~RAD_OpticalDepthZ () {return;}
  Real getOpticalDepth (Real z1, Real z2);
  double integrand (double z);
private:
  RAD_OpticalDepth* itsRADOD;
  // to prevent copying and assignment:
  RAD_OpticalDepthZ (const RAD_OpticalDepthZ & A);
  RAD_OpticalDepthZ operator = (const RAD_OpticalDepthZ& A);
};

class RAD_OpticalDepthU : public Integral {
public:
  RAD_OpticalDepthU (RAD_OpticalDepth* A);
  ~RAD_OpticalDepthU () {return;}
  Real getOpticalDepth (Real z1, Real z2);
  Real getOpticalDepth (Real z);
  double integrand (double z);
private:
  bool isPositiveMu;
  RAD_OpticalDepth* itsRADOD;
  // to prevent copying and assignment:
  RAD_OpticalDepthU (const RAD_OpticalDepthU & A);
  RAD_OpticalDepthU operator = (const RAD_OpticalDepthU& A);
};

class RAD_OpticalDepth {
 public:
  RAD_OpticalDepth (Velocity* V, Real DeltaE, Real Gamma, Real Tau0,
		    Real Vinfty); // constructor
  //~RAD_OpticalDepth (); // destructor
  Real getOpticalDepth (Real p, Real z);
  Real getPhi (Real wz);
  Real getP () {return itsP;}
  Real getW (Real u) {return itsVelocity->getVelocity (u);}
  void initialize (Real p, Real z); // setup to manually call integrand
  // initialize is also called by getOpticalDepth
  double getZIntegrand (double z);
 private:
  static const Real LARGE_OPTICAL_DEPTH;
  static const Real MU0;
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
  bool isTransparent;
  Velocity* itsVelocity;
  RAD_OpticalDepthZ* itsRADODZ;
  RAD_OpticalDepthU* itsRADODU;
  void checkInput ();
  void allocateRAD_OpticalDepthZU ();
  void freeRAD_OpticalDepthZU ();
  // To prevent copying and assignment:
  //RAD_OpticalDepth (RAD_OpticalDepth const & T);
  //RAD_OpticalDepth operator = (RAD_OpticalDepth const & T);

};

# endif//RAD_OPTICAL_DEPTH_H
