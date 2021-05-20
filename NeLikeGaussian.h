/***************************************************************************
    NeLikeGaussian.h -  a class to get gaussian lines for all three 
                        2p-3s lines in a Ne-like ion 

                             -------------------
    begin				: May 2021
    copyright			: (C) 2021 by Maurice Leutenegger
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

#ifndef NE_LIKE_GAUSSIAN_H
#define NE_LIKE_GAUSSIAN_H

#include "xsTypes.h"

class NeLikeGaussian
{
 public:
  NeLikeGaussian (const RealArray& energy, const RealArray& parameter);
  /* initialization assumes inclusion of calibration parameter 
     (delta lambda mA) and 3F shift paramteter (mA) */
  ~NeLikeGaussian ();
  void getFlux (RealArray& flux);
 private:
  const static size_t Nlines = 3; 
  const RealArray itsEnergyArray;
  Real its3F; // ratio of 3F / (3F+3G+M2)
  Real itsM2_3G; // ratio of M2/3G
  Real itsSigma;
  RealArray itsEnergy; 
  void extractParameters (const RealArray& parameter);
};

#endif//HE_LIKE_GAUSSIAN_H
