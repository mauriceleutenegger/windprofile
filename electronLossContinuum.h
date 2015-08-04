/***************************************************************************
    electronLossContinuum.h   - compute an electron loss continuum model 
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

#ifndef ELECTRONLOSSCONTINUUM_H
#define ELECTRONLOSSCONTINUUM_H

#include "xsTypes.h"
#include "mal_integration.h"

class electronLossContinuum : public Integral
{
 public :
  electronLossContinuum (Real sigma, Real ETau, Real elc_tail_ratio, Real E0);
  ~electronLossContinuum ();
  Real getELC (Real deltaE); /* deltaE = E0 - E */
  /* Note that for deltaE < 0, you can still get some flux from the
   Gaussian wings of the convolution. */
  double integrand (double epsilon);
 private:
  void checkInput ();
  Real itsSigma;
  Real itsETau;
  Real itsELCTailRatio;
  Real itsE0;
  Real itsDeltaE; 
};
#endif /* ELECTRONLOSSCONTINUUM_H */
