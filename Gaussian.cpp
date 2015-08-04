/***************************************************************************
    Gaussian.cpp - computes the a Gaussian profile for XSPEC.

                             -------------------
    begin				: August 2007
    copyright			: (C) 2007 by Maurice Leutenegger
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

#include "Gaussian.h"
#include "Utilities.h"
#include <iostream>
#include <math.h>
/*
  ML - Oct 16 2013
  The implementation of the calls to erfc used to be through gsl - 
  specifically <gsl/gsl_sf_erf.h> -
  but there is a perfectly valid implementation in <math.h>. 
  As far as I know there is no reason to introduce a dependency on gsl,
  so I have cut it out.
*/


using namespace std;

const Real Gaussian::itsNumberOfSigmas = 6.0;

Gaussian::Gaussian (const RealArray& energy, Real LineEnergyKEV, Real SigmaC) 
  : itsEnergyArray (energy), itsLineEnergy (LineEnergyKEV), itsSigma (SigmaC), 
    itsEnergySize (energy.size()), itsFluxSize (itsEnergySize - 1), 
    itsCentralBin (itsFluxSize / 2), itsMinimumBin (0), 
    itsMaximumBin (itsFluxSize), isBackwards (true), isOutsideRange (true), 
    isUnresolved (true)
{
  checkInput ();
  return;
}

Gaussian::~Gaussian ()
{
  return;
}

void Gaussian::checkInput ()
{
  /* 1. Check response direction. This code assumes an ascending response. */
  if (itsEnergyArray[1] > itsEnergyArray[0])
    isBackwards = false;
  else isBackwards = true;
  /* 2. Set and check range. */
  itsCentralBin = BinarySearch (itsEnergyArray, itsLineEnergy);
  Real MinimumEnergy = itsLineEnergy * (1. - itsNumberOfSigmas * itsSigma);
  Real MaximumEnergy = itsLineEnergy * (1. + itsNumberOfSigmas * itsSigma);
  isOutsideRange = ((itsLineEnergy < itsEnergyArray[0]) || (itsEnergyArray[itsEnergySize - 1] < itsLineEnergy));
  if (!isOutsideRange) {
    itsMinimumBin = BinarySearch (itsEnergyArray, MinimumEnergy); 
    itsMaximumBin = BinarySearch (itsEnergyArray, MaximumEnergy);
  }
  /* 3. Check that sigma > 0. */
  isUnresolved = true;
  if (itsSigma > 0.) isUnresolved = false;
  /* Check if the line is effectively unresolved. */
  if ((itsMaximumBin - itsMinimumBin) < 3) isUnresolved = true;
  return;
}

void Gaussian::setParameters (Real LineEnergyKEV, Real SigmaC)
{
  itsLineEnergy = LineEnergyKEV;
  itsSigma = SigmaC;
  checkInput ();
  return;
}

void Gaussian::getFlux (RealArray& flux)
{
  /* 1. Initialize flux. */
  flux.resize(itsFluxSize);
  /* 2. If response is backwards, bail. */
  if (isBackwards) { 
    cout << "Gaussian::getFlux: Backwards response not supported.\n";
    return;
  }
  /* 3. If line is outside range, bail */
  if (isOutsideRange) { 
    cout << "Gaussian::getFlux: rest energy (" << itsLineEnergy << 
      " keV) is outside response range.\n";
    return;
  }
  /* 4. If width is zero, put flux in central bin. */
  if (isUnresolved) { 
    Real binCenterEnergy = (itsEnergyArray[itsCentralBin+1] + 
			    itsEnergyArray[itsCentralBin]) / 2.;
    Real binWidth = itsEnergyArray[itsCentralBin+1] - 
      itsEnergyArray[itsCentralBin];
    Real centeringFraction = (itsLineEnergy - binCenterEnergy) / binWidth;
    flux[itsCentralBin] = 1. - fabs (centeringFraction);
    if (centeringFraction > 0.) {
      if (itsCentralBin < (itsFluxSize - 1)) {
	flux[itsCentralBin + 1] = centeringFraction;
      }
    } else { /* centeringFraction < 0. */
      if (itsCentralBin > 0) {
	flux[itsCentralBin - 1] = fabs (centeringFraction);
      }
    }
    return;
  }
  /* 5. Calculate Gaussian:
     The Gaussian is split into two parts so that each half can be computed 
     in such a way as to minimize rounding errors.
  */
  RealArray GaussianArgument (itsEnergySize);
  GaussianArgument = 
    (itsEnergyArray / itsLineEnergy - 1.) / (itsSigma * M_SQRT2);

  Real previous = erfc (fabs(GaussianArgument[itsMinimumBin]));
  for (size_t i = itsMinimumBin; i < itsCentralBin; i++) {
    Real current = erfc (fabs(GaussianArgument[i+1]));
    flux[i] = 0.5 * fabs (current - previous); 
    previous = current;
  }

  previous = erfc (GaussianArgument[itsMaximumBin]);
  for (size_t i = itsMaximumBin; i > itsCentralBin; i--) {
    Real current = erfc (GaussianArgument [i-1]);
    flux[i-1] = 0.5 * fabs (current - previous); 
    previous = current;
  }
  return;
}

