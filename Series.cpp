/***************************************************************************
    Series.h   - Base type for iterative calculation of a taylor series.

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

#include "Series.h"
#include "Utilities.h"
#include <iostream>

using namespace std;

Series::Series (double tolerance) : itsTolerance (tolerance), n (0)
{
  checkInput ();
  return;
}

Series::~Series ()
{
  return;
}

void Series::checkInput ()
{
  if (compare (itsTolerance, 0.1) == 1) itsTolerance = 0.1;
  if (compare (itsTolerance, 0.) == -1) itsTolerance = 0.1;
  return;
}

double Series::sumSeries ()
{
  double sum = getTerm ();
  double term = 0.;
  double error = 0.;
  n = 1;
  do {
    iterate ();
    term = getTerm ();
    sum += term;
    error = fabs (term / sum);
    ++n;
    if (n > MAXIMUM_ITERATIONS) {
      cout << "Series: convergence not achieved in " << MAXIMUM_ITERATIONS << 
	"; aborting.\n";
      break;
    }
  } while (compare (error, itsTolerance) == 1);
  return sum;
}

