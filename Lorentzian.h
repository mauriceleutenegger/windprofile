/***************************************************************************
    Lorentzian.h - computes a Lorentzian profile for XSPEC.

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

#ifndef LORENTZIAN_H
#define LORENTZIAN_H

#include "xsTypes.h"

class Lorentzian {
 public:
  Lorentzian (const RealArray& Energy, Real LineEnergyKEV, Real GammaC);
  ~Lorentzian ();
  void setParameters (Real LineEnergyKEV, Real GammaC);
  void getFlux (RealArray& flux);
  bool checkOutOfBounds () {return isOutsideRange;}
 private:
  static const Real itsNumberOfGammas;
  const RealArray& itsEnergyArray;
  Real itsLineEnergy;
  Real itsGamma;
  size_t itsEnergySize;
  size_t itsFluxSize;
  size_t itsCentralBin;
  size_t itsMinimumBin;
  size_t itsMaximumBin;
  bool isBackwards; /* True if energy decreases with increasing index. */
  bool isOutsideRange; /* True if line center is outside response range. */
  bool isUnresolved; /* True if gamma = 0. */
  void checkInput ();
  Real getLorentzianArgument (int n);
};

#endif//LORENTZIAN_H
