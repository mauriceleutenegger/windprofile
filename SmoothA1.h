/***************************************************************************
    SmoothA1.h   - Calculates part of the analytic expression for
                   smooth optical depth using a series expansion.

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

#ifndef SMOOTH_A1
#define SMOOTH_A1

#include "Series.h"

class SmoothA1 : public Series {
 public:
  SmoothA1 (double p = 1., double z = 1.);
 private:
  double getTerm ();
  void iterate ();
  void initialize ();
  void checkInput ();
  double p;
  double z;
  double mu2;
  double muTerm;
  double z2;
  double zTerm;
  double zStar2;
  double zStarTerm;
};

#endif//SMOOTH_A1
