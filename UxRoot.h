/***************************************************************************
    UxRoot.h   - Finds the u coordinate at which p = 1 + epsilon for a given
                 value of x. This allows the avoidance of the occulted part
                 of the Lx integral.

                             -------------------
    begin				: March 2007
    copyright			: (C) 2007 by Maurice Leutenegger
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

#ifndef UXROOT_H
#define UXROOT_H

#include "mal_RootFinderNewton.h"
#include "Utilities.h"
#include "xsTypes.h"

class UxRoot : public RootFinderNewton {
 public:
  UxRoot (Velocity* V, Real x);
  double f (double u);
  double df (double u);
  void fdf (double u, double& y, double& dy);
  void setX (Real x);
 private:
  Velocity* itsVelocity;
  Real itsX2;
};


#endif//UXROOT_H
