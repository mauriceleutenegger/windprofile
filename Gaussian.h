/***************************************************************************
    Gaussian.h - computes the a Gaussian profile for XSPEC.

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

#ifndef GAUSSIAN_H
#define GAUSSIAN_H

#include "xsTypes.h"

class Gaussian {
 public:
  Gaussian (const RealArray& Energy, Real LineEnergyKEV, Real SigmaC);
  ~Gaussian ();
  void setParameters (Real LineEnergyKEV, Real SigmaC);
  void getFlux (RealArray& flux);
  bool checkOutOfBounds () {return isOutsideRange;}
 private:
  static const Real itsNumberOfSigmas;
  const RealArray& itsEnergyArray;
  Real itsLineEnergy;
  Real itsSigma;
  size_t itsEnergySize;
  size_t itsFluxSize;
  size_t itsCentralBin;
  size_t itsMinimumBin;
  size_t itsMaximumBin;
  bool isBackwards; /* True if energy decreases with increasing index. */
  bool isOutsideRange; /* True if line center is outside response range. */
  bool isUnresolved; /* True if sigma = 0. */
  void checkInput ();
  Real getGaussianArgument (int n);
};

#endif//GAUSSIAN_H
