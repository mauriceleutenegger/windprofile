/***************************************************************************
    HeLikeGaussian.h -  a class to get gaussian lines for all lines in
                        a He-like triplet (including both intercombo lines).

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

#ifndef HE_LIKE_GAUSSIAN_H
#define HE_LIKE_GAUSSIAN_H

#include "xsTypes.h"

class HeLikeGaussian
{
 public:
  HeLikeGaussian (const RealArray& energy, const RealArray& parameter);
  /* initialization assumes inclusion of calibration parameter 
     (delta lambda mA) */
  ~HeLikeGaussian ();
  void getFlux (RealArray& flux);
 private:
  const static size_t Nlines = 4;
  const RealArray itsEnergyArray;
  Real itsR;
  Real itsG;
  Real itsSigma;
  RealArray itsEnergy;
  Real itsXFraction;
  bool isXFinite;
  void extractParameters (const RealArray& parameter);
};

#endif//HE_LIKE_GAUSSIAN_H
