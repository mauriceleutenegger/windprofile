/***************************************************************************
    IsotropicSeries.cpp   - Calculates part of the analytic expression for 
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

#include "IsotropicSeries.h"
#include "Utilities.h"
#include <iostream>

using namespace std;

IsotropicSeries::IsotropicSeries (double x) // x = zh / z
  : Series (), x2 (x * x), xterm (1.)
{
  checkInput ();
  return;
}

void IsotropicSeries::checkInput ()
{
  if (compare (x2, 1.) != -1) {
    cerr << "IsotropicSeries: x too large. x2 = " << x2 << "\n";
    x2 = 0.;
  }
  return;
}

double IsotropicSeries::getTerm ()
{
  return (xterm / (2. * getN () + 1.));
}

void IsotropicSeries::iterate ()
{
  xterm *= -1. * x2;
  return;
}
