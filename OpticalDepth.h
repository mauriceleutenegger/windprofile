/***************************************************************************
    OpticalDepth.h   - Calculates the optical depth to X-rays for a stellar
                       wind along a ray p from point z.

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

#ifndef OPTICAL_DEPTH_H
#define OPTICAL_DEPTH_H

#include <stdbool.h>
#include "xsTypes.h"
#include "Utilities.h"
#include "Porosity.h"
#include "AnalyticOpticalDepth.h"
#include "NumericalOpticalDepth.h"

/***********************************************************************\

Upon instantiation, this class allocates memory for either an 
AnalyticOpticalDepth, or for a NumericalOpticalDepth along with a Porosity
and a Velocity. The destructor handles the deletion.

Note that although taustar and h may be changed after instantiation,
so that one may grid over parameters, beta must remain constant,
and the booleans are fixed. Thus, to study different optical depth 
formulations, one must instantiate a new OpticalDepth.

\***********************************************************************/
class OpticalDepth {
 public:
  /*  OpticalDepth 
    (Real TauStar = 0., Real h = 0., Real beta = 1., bool numerical = false,
    bool anisotropic = false, bool rosseland = false, bool expansion = false);*/
  OpticalDepth 
    (Real TauStar = 0., Real h = 0., Real beta = 1., bool numerical = false,
     bool anisotropic = false, bool prolate = false,
     bool rosseland = false, bool expansion = false, 
     bool HeII = false);
  ~OpticalDepth ();
  Real getOpticalDepth (Real p, Real z);
  void setParameters (Real TauStar, Real h);
 private:
  static const Real MINIMUM_VELOCITY; // scaled velocity at R*
  bool isNumerical;
  bool isHeII;
  Velocity* itsVelocity;
  Porosity* itsPorosity;
  AnalyticOpticalDepth* itsAnalyticOpticalDepth;
  NumericalOpticalDepth* itsNumericalOpticalDepth;
  // allocators and deallocators
  void allocateVelocity (Real beta);
  void allocatePorosity 
    (Real TauStar, Real h, bool anisotropic, bool prolate, bool rosseland);
  void allocateNumericalOpticalDepth (Real TauStar);
  void allocateAnalyticOpticalDepth 
    (Real TauStar, Real h, bool anisotropic, bool expansion);
  void freeClasses ();
  // To prevent copying and assignment:
  OpticalDepth (OpticalDepth const & T);
  OpticalDepth operator = (OpticalDepth const & T);
};

#endif//OPTICAL_DEPTH_H
