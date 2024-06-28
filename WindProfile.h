/***************************************************************************
    WindProfile.h   - This class computes the flux array for the XSPEC models
                      by integrating Lx dx with FluxIntegral. It uses
                      WindParameter to set the model parameters.

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

#ifndef WIND_PROFILE_H
#define WIND_PROFILE_H

#include "xsTypes.h"
#include "Utilities.h"
#include "Porosity.h"
#include "OpticalDepth.h"
#include "NumericalOpticalDepth.h"
#include "HeLikeRatio.h"
#include "ResonanceScattering.h"
#include "Lx.h"
#include "WindParameter.h"
#include "FluxIntegral.h"

class WindProfile
{
 public:
  WindProfile (const RealArray& energy, const RealArray& parameter, 
	       ModelType type = general);
  ~WindProfile ();
  void getModelFlux (RealArray& flux);
 private:
  const RealArray& itsEnergyArray;
  size_t itsEnergySize;
  size_t itsFluxSize;
  RealArray x;
  ModelType itsModelType;
  WindParameter* itsWindParameter;
  FluxIntegral* itsFluxIntegral;
  Lx* itsLx;
  Velocity* itsVelocity;
  HeLikeRatio* itsHeLikeRatio;
  ResonanceScattering* itsResonanceScattering;
  OpticalDepth* itsOpticalDepth;
  OpticalDepth* itsOpticalDepthHeII; // NEW
  /* I am pretty sure that the NumericalOpticalDepth declared here
   never gets used. */
  //  NumericalOpticalDepth* itsNumericalOpticalDepth;
  Porosity* itsPorosity;
  RAD_OpticalDepth* itsRAD_OpticalDepth;
  Real itsTotal;
  bool isFinite;
  bool isNumerical;
  bool isHeII;
  void allocateWindParameter (const RealArray& parameter);
  void freeWindParameter ();
  void allocateClasses ();
  void freeClasses ();
  void getOneFlux (RealArray& flux, HeLikeType type = wResonance);
  void renormalize (RealArray& flux);
  void TransmissionRatio (const RealArray& x);
  void FToIRatio (const RealArray& fFlux, const RealArray& iFlux);
};

#endif//WIND_PROFILE_H
