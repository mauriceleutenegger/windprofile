/***************************************************************************
    IsotropicSeries.h   - Calculates part of the analytic expression for 
                          optical depth with isotropic clumps using a 
                          series expansion.

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

#ifndef ISOTROPIC_SERIES
#define ISOTROPIC_SERIES

#include "Series.h"

class IsotropicSeries : public Series {
 public:
  IsotropicSeries (double x = 0.);
 private:
  double x2;
  double xterm;
  double getTerm ();
  void iterate ();
  void checkInput ();
};

#endif//ISOTROPIC_SERIES
