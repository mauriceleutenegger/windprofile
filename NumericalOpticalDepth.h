/***************************************************************************
    NumericalOpticalDepth.h   - Computes the optical depth numerically
                                by breaking the integral into pieces and
                                integrating the pieces along du or dz.

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

#ifndef NUMERICAL_OPTICAL_DEPTH_H
#define NUMERICAL_OPTICAL_DEPTH_H

#include <stdbool.h>
#include <gsl/gsl_math.h>
#include "xsTypes.h"
#include "Utilities.h"
#include "Porosity.h"
#include "mal_integration.h"

class NumericalOpticalDepth;

class NumericalOpticalDepthZ : public Integral {
 public:
  NumericalOpticalDepthZ (NumericalOpticalDepth* A);
  ~NumericalOpticalDepthZ () {return;}
  Real getOpticalDepth (Real p, Real z1, Real z2);
  Real getOpticalDepth (Real p, Real z);
  Real getOpticalDepth (Real p);
  double integrand (double z);
 private:
  Real itsP;
  NumericalOpticalDepth* itsNumericalOpticalDepth;
  // to prevent copying and assignment:
  NumericalOpticalDepthZ (const NumericalOpticalDepthZ& A);
  NumericalOpticalDepthZ operator = (const NumericalOpticalDepthZ& A);
};

class NumericalOpticalDepthU : public Integral {
 public:
  NumericalOpticalDepthU (NumericalOpticalDepth* A);
  ~NumericalOpticalDepthU () {return;}
  Real getOpticalDepth (Real p, Real z);
  double integrand (double u);
 private:
  Real itsP;
  NumericalOpticalDepth* itsNumericalOpticalDepth;
  // to prevent copying and assignment:
  NumericalOpticalDepthU (const NumericalOpticalDepthU& A);
  NumericalOpticalDepthU operator = (const NumericalOpticalDepthU& A);
};

class NumericalOpticalDepth 
{
 public:
  /*  NumericalOpticalDepth
      (Real TauStar, Porosity* P, Velocity* V);*/
  NumericalOpticalDepth
    (Real TauStar, bool HeII, Porosity* P, Velocity* V);
  ~NumericalOpticalDepth ();
  void setParameters (Real TauStar, Porosity* P, Velocity* V);
  void setTauStar (Real TauStar);
  Real getOpticalDepth (Real p, Real z);
  friend double NumericalOpticalDepthZ::integrand (double z);
  friend double NumericalOpticalDepthU::integrand (double u);
 private:
  static const Real LARGE_OPTICAL_DEPTH;
  static const Real MU0;
  static const Real LARGE_P;
  Real itsP; // impact parameter p
  Real itsTauStar; // characteristic optical depth
  bool isTransparent;
  bool isPorous;
  Porosity* itsPorosity; 
  Velocity* itsVelocity;
  NumericalOpticalDepthZ* itsNumericalOpticalDepthZ;
  NumericalOpticalDepthU* itsNumericalOpticalDepthU;
  void checkInput ();
  void allocateNumericalOpticalDepthZU ();
  void freeNumericalOpticalDepthZU ();
  // for He II optical depth
  bool isHeII;
  // The Butterworth filter parameters are hardcoded for now.
  int nButterworth; 
  Real uButterworth; 
  Real HeIIFilter (Real u); 
  // To prevent copying or assignment:
  NumericalOpticalDepth (const NumericalOpticalDepth& Tau);
  NumericalOpticalDepth operator = (const NumericalOpticalDepth& Tau);
};

#endif//NUMERICAL_OPTICAL_DEPTH_H
