/***************************************************************************
    calorimeterLSF.cpp   - XSPEC models to compute calorimeter LSFs.
                     

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

#ifndef CALORIMETERLSF_H
#define CALORIMETERLSF_H

#include "xsTypes.h"
#include "mal_integration.h"
#include "electronLossContinuum.h"

class calorimeterLSF : public Integral
{
 public:
  calorimeterLSF (Real sigma, Real eTau, Real elc_tail_ratio, Real E0);
  ~calorimeterLSF ();
  Real getLSF (Real deltaE1, Real deltaE2);
  double integrand (double deltaE);
 private:
  void allocateClasses ();
  void freeClasses ();
  void checkInput ();
  Real itsSigma; // Gaussian width of core LSF
  Real itsETau; // Characteristic energy scale of electron loss continuuum
  Real itsELCTailRatio;
  Real itsE0;
  electronLossContinuum* itsELC;
};

#endif /* CALORIMETERLSF_H */
