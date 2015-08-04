/***************************************************************************
    Lorentzian.cpp - computes a Lorentzian profile for XSPEC.

                             -------------------
    begin				: October 2013
    copyright			: (C) 2013 by Maurice Leutenegger
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

#include "Lorentzian.h"
#include "Utilities.h"
#include <iostream>
#include <math.h>

using namespace std;

const Real Lorentzian::itsNumberOfGammas = 5.0;

Lorentzian::Lorentzian (const RealArray& energy, Real LineEnergyKEV, Real GammaC) 
  : itsEnergyArray (energy), itsLineEnergy (LineEnergyKEV), itsGamma (GammaC), 
    itsEnergySize (energy.size()), itsFluxSize (itsEnergySize - 1), 
    itsCentralBin (itsFluxSize / 2), itsMinimumBin (0), 
    itsMaximumBin (itsFluxSize), isBackwards (true), isOutsideRange (true), 
    isUnresolved (true)
{
  checkInput ();
  return;
}

Lorentzian::~Lorentzian ()
{
  return;
}

void Lorentzian::checkInput ()
{
  /* 1. Check response direction. */
  if (itsEnergyArray[1] > itsEnergyArray[0])
    isBackwards = false;
  else isBackwards = true;
  /* 2. Set and check range. */
  itsCentralBin = BinarySearch (itsEnergyArray, itsLineEnergy);
  isOutsideRange = ((itsLineEnergy < itsEnergyArray[0]) || (itsEnergyArray[itsEnergySize - 1] < itsLineEnergy));
  /* 3. Check that gamma > 0. */
  isUnresolved = true;
  if (itsGamma > 0.) isUnresolved = false;
  return;
  /* Note that minimum and maximum bin are currently unused. I am not 
     bothering to check for aliasing problems in the core of the profile 
     because I am assuming that anytime you would bother to use a Lorentzian 
     instead of a Gaussian, it's because the Lorentzian is broad enough 
     that you care about the wings. */
}

void Lorentzian::setParameters (Real LineEnergyKEV, Real GammaC)
{
  itsLineEnergy = LineEnergyKEV;
  itsGamma = GammaC;
  checkInput ();
  return;
}

void Lorentzian::getFlux (RealArray& flux)
{
  /* 1. Initialize flux. */
  flux.resize(itsFluxSize, 0.);
  /* 2. If response backwards, bail. */
  if (isBackwards) { 
    cout << "Lorentzian::getFlux: Backwards response not supported.\n";
    return;
  }
  /* 3. If line outside range, bail. */
  if (isOutsideRange) { 
    cout << "Lorentzian::getFlux: rest energy (" << itsLineEnergy << 
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
  /* 5. Calculate Lorentzian. */

  /*
    The Lorentzian is split into two parts so that each half can be computed 
    in such a way as to minimize rounding errors.
  */
  RealArray LorentzianArgument (itsEnergySize);
  LorentzianArgument = 
    (itsEnergyArray / itsLineEnergy - 1.) / itsGamma;

  /* 
     Note that the loops go out to the end of the response matrix, 
     unlike the Gaussian,
     since the Lorentzian wings are broad, and the computation 
     of the profile is not expensive. However, at some point
     if the Lorentzian is sufficiently broad, the idea of a
     "non-relativistic" velocity shift is not valid.
  */
  Real previous = atan (fabs(LorentzianArgument[0]));
  for (size_t i = 0; i < itsCentralBin; i++) {
    Real current = atan (fabs(LorentzianArgument[i+1]));
    flux[i] = fabs (current - previous) / M_PI; 
    previous = current;
  }

  previous = atan (LorentzianArgument[itsEnergySize - 1]);
  for (size_t i = itsFluxSize; i > itsCentralBin; i--) {
    Real current = atan (LorentzianArgument [i-1]);
    flux[i-1] = fabs (current - previous) / M_PI; 
    previous = current;
  }

  return;
}

