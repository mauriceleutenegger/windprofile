/***************************************************************************
    AnalyticOpticalDepth.h   - Computes wind optical depth analytically.

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

#ifndef ANALYTIC_OPTICAL_DEPTH_H
#define ANALYTIC_OPTICAL_DEPTH_H

#include <stdbool.h>
#include <float.h>
#include "xsTypes.h"
#include "SmoothA1.h"
#include "IsotropicSeries.h"

/* Analytic optical depth, beta = 1. Includes a few porosity options:
   Smooth, isotropic expansion, isotropic stretch, anisotropic stretch. The 
   isotropic models use a Rosseland form for the bridging law, while the 
   anisotropic model has a step function. */
class AnalyticOpticalDepth 
{
 public:
  AnalyticOpticalDepth 
  (Real TauStar = 0., Real h = 0., bool anisotropic = false, 
   bool expansion = false);
  Real getOpticalDepth (Real p, Real z);
  void setParameters (Real TauStar, Real h);
 private:
  static const Real LARGE_OPTICAL_DEPTH;
  // model parameters
  Real itsTauStar;
  Real itsTauClump0; // = TauStar * h
  bool isTransparent;
  bool isPorous;
  bool isAnisotropic;
  bool isExpansion;
  // variables
  Real mu; 
  Real zstar;
  //
  Real Smooth (Real p, Real z);
  Real SmoothG1 (Real p, Real z);
  Real SmoothL1 (Real p, Real z);
  Real Expansion (Real p, Real z);
  Real ExpansionG1 (Real p, Real z);
  Real Expansion1 (Real p, Real z);
  Real ExpansionL1 (Real p, Real z);
  Real Isotropic (Real p, Real z);
  Real Isotropic1 (Real p, Real z);
  Real Anisotropic (Real p, Real z);
  Real Anisotropic1 (Real ra, Real rb);
  void checkInput ();
};

#endif//ANALYTIC_OPTICAL_DEPTH_H
