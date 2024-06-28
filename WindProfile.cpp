/***************************************************************************
    WindProfile.cpp   - This class computes the flux array for the XSPEC models
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

#include "WindProfile.h"
#include <iostream>

using namespace std;

WindProfile::WindProfile 
(const RealArray& energy, const RealArray& parameter, ModelType type)
  : itsEnergyArray (energy),  itsEnergySize (itsEnergyArray.size ()), 
    itsFluxSize (itsEnergySize - 1), x (RealArray (itsEnergySize)), 
    itsModelType (type), itsWindParameter (NULL), itsFluxIntegral (NULL),
    itsLx (NULL), itsVelocity (NULL), itsHeLikeRatio (NULL),
    itsResonanceScattering (NULL), itsOpticalDepth (NULL),
    //    itsNumericalOpticalDepth (NULL),
    itsPorosity (NULL),
    itsRAD_OpticalDepth (NULL),
    itsTotal (0.), 
    isFinite (false), isNumerical (false)
{
  allocateWindParameter (parameter);
  allocateClasses ();
  return;
}

WindProfile::~WindProfile ()
{
  freeWindParameter ();
  freeClasses ();
}

void WindProfile::allocateWindParameter (const RealArray& parameter)
{
  itsWindParameter = new WindParameter (parameter, itsModelType);
  return;
}

void WindProfile::freeWindParameter ()
{
  delete itsWindParameter;
  itsWindParameter = NULL;
  return;
}

void WindProfile::allocateClasses ()
{
  isNumerical = itsWindParameter->getNumerical ();
  isHeII = itsWindParameter->getHeII ();
  itsWindParameter->initializeVelocity (itsVelocity);
  if (itsModelType == helike) {
    itsWindParameter->initializeHeLikeRatio (itsHeLikeRatio);
  } else itsHeLikeRatio = 0;
  // HeLikeRatio will only be utilized if it's not a pointer to NULL.
  itsWindParameter->initializeResonanceScattering 
    (itsResonanceScattering, itsVelocity);
  itsWindParameter->initializeOpticalDepth (itsOpticalDepth, 
					    itsOpticalDepthHeII);
  if (itsModelType == rad) {
    itsWindParameter->initializeRAD_OpticalDepth (itsRAD_OpticalDepth,
						 itsVelocity);
  }
  if (isHeII) {
    itsWindParameter->initializeLx (itsLx, itsVelocity, itsHeLikeRatio, 
				    itsResonanceScattering, itsOpticalDepth,
				    itsOpticalDepthHeII);
  } else if (itsModelType == rad) {
    itsWindParameter->initializeLx (itsLx, itsVelocity, itsHeLikeRatio, 
				    itsResonanceScattering, itsOpticalDepth,
				    itsRAD_OpticalDepth);
  } else {
    itsWindParameter->initializeLx (itsLx, itsVelocity, itsHeLikeRatio, 
				    itsResonanceScattering, itsOpticalDepth);
  }
  itsFluxIntegral = new FluxIntegral (itsLx);
  return;
}

// should be OK to delete NULL on optical depth
void WindProfile::freeClasses ()
{
  delete itsWindParameter;
  itsWindParameter = NULL;
  delete itsFluxIntegral;
  itsFluxIntegral = NULL;
  delete itsLx;
  itsLx = NULL;
  delete itsOpticalDepth;
  itsOpticalDepth = NULL;
  if (isHeII) {
    delete itsOpticalDepthHeII;
    itsOpticalDepthHeII = NULL;
  }
  //  delete itsNumericalOpticalDepth;
  //  itsNumericalOpticalDepth = NULL;
  delete itsPorosity;
  itsPorosity = NULL;
  delete itsResonanceScattering;
  itsResonanceScattering = NULL;
  delete itsHeLikeRatio;
  itsHeLikeRatio = NULL;
  delete itsVelocity;
  itsVelocity = NULL;
  return;
}

void WindProfile::getModelFlux (RealArray& flux) 
{
  flux.resize (itsFluxSize);
  Real RADeff = 1.;
  if (itsModelType == helike) {
    RealArray rFlux (itsFluxSize);
    RealArray iFlux (itsFluxSize);
    RealArray fFlux (itsFluxSize);
    getOneFlux (rFlux, wResonance);
    getOneFlux (iFlux, yIntercombination);
    getOneFlux (fFlux, zForbidden);
    if (itsWindParameter->getVerbosity ()) {
      FToIRatio (fFlux, iFlux);
    }
    Real G = itsWindParameter->getG ();
    flux = (rFlux + G * (iFlux + fFlux)) / (1. + G);
  } else if (itsModelType == rad) {
    RealArray noRADflux (itsFluxSize);
    getOneFlux (flux);
    itsLx->setRADTransparent ();
    getOneFlux (noRADflux);
    itsLx->notRADTransparent ();
    Real noRADflux_sum = noRADflux.sum ();
    if (compare (noRADflux_sum, 0.) == 1) {
      RADeff = flux.sum () / noRADflux.sum ();
    }
  } else {
    getOneFlux (flux);
  }
  renormalize (flux);
  if (itsModelType == rad) {
    flux *= RADeff;
    cout << "RAD transmitted fraction: " << RADeff << endl;
  }
  if (isFinite && itsWindParameter->getVerbosity () && 
      (itsModelType != helike)) {
    TransmissionRatio (x);
  }
  return;
}

void WindProfile::getOneFlux (RealArray& flux, HeLikeType type)
{
  if (itsModelType == helike) {
    itsWindParameter->setX (itsEnergyArray, x, type);
    itsLx->setHeLikeType (type);
  } else {
    itsWindParameter->setX (itsEnergyArray, x);
  }
  for (size_t i = 0; i < (itsFluxSize); i++) {
    flux[i] = itsFluxIntegral->getFlux (x[i], x[i+1]);
  } 
  return;
}

void WindProfile::renormalize (RealArray& flux)
{
  itsTotal = flux.sum ();
  if (compare (itsTotal, 0.) == 1) {
    flux /= itsTotal;
    isFinite = true;
    return;
  } else {
    cerr << "Can't renormalize; total flux is zero.\n";
    isFinite = false;
    return;
  }
}

void WindProfile::TransmissionRatio (const RealArray& x)
{
  itsLx->setTransparent ();
  Real UnabsorbedTotal = 0.;
  for (size_t i = 0; i < itsFluxSize; i++) {
    UnabsorbedTotal += itsFluxIntegral->getFlux (x[i], x[i+1]);
  }
  if (compare (UnabsorbedTotal, 0.) == 1) {
    Real ratio = itsTotal / UnabsorbedTotal;
    cout << "(observed photons) / (emitted photons) = " << ratio << "\n";
  }
  return;
}

void WindProfile::FToIRatio (const RealArray& fFlux, const RealArray& iFlux)
{
  if (compare (iFlux.sum (), 0.) != 1) {
    cout << "Can't compute f/i ratio for iFlux = 0.\n";
    return;
  }
  Real ratio = fFlux.sum () / iFlux.sum ();
  cout << "R = f / i = " << ratio << "\n";
  return;
}
