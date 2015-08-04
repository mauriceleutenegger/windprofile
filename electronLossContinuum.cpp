/***************************************************************************
    electronLossContinuum.cpp   - compute an electron loss continuum model 
                                convolved with a gaussian.
                     

                             -------------------
    begin				: November 2013
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

#include "electronLossContinuum.h"
#include <iostream>
#include "Utilities.h"

electronLossContinuum::electronLossContinuum 
(Real sigma, Real ETau, Real elc_tail_ratio, Real E0) 
  : itsSigma (sigma), itsETau (ETau), itsELCTailRatio (elc_tail_ratio), itsE0 (E0), itsDeltaE (0.)
{
  checkInput ();
  return;
}

electronLossContinuum::~electronLossContinuum ()
{
  return;
}

void electronLossContinuum::checkInput ()
{
  if (itsSigma < 0.) {
    itsSigma = 0.;
    cerr << "electronLossContiuum: sigma < 0, setting to 0\n";
  }
  return;
}

/* epsilon is the energy in the convolution integral */
double electronLossContinuum::integrand (double epsilon) 
{
  double tail = 0.;
  if (compare (itsDeltaE, epsilon) == 1) {
    tail = exp (-1. * (itsDeltaE - epsilon) / itsETau) / itsETau;
    //    tail = exp (-1. * (itsDeltaE - epsilon) / itsETau);
    /* Add in a constant floor. */
    /* Needs to be in units of 1/E also! */
    tail += itsELCTailRatio * itsE0;

    // binFraction is the width of the energy bin divided by the line energy
    // need to know energy bin width
    //    tail += itsNFlat;
    //    tail /= itsETau;
    // need to redo the normalization...
  }
  else
    return 0.;
  double gnorm = 1. / (itsSigma * sqrt (2 * M_PI));
  double arg = epsilon / itsSigma;
  double gaussian = gnorm * exp (-0.5 * arg * arg);
  return tail * gaussian;
}

Real electronLossContinuum::getELC (Real deltaE)
{
  itsDeltaE = deltaE;
  int Nsigma = 10;
  Real bound = Nsigma * itsSigma;
  return qags (-1 * bound, bound);
  //  return qag (-1 * bound, bound);
  //  return qagi ();
}
