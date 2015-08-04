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

#include "xsTypes.h"
#include "calorimeterLSF.h"
#include "Utilities.h"
#include <iostream>

calorimeterLSF::calorimeterLSF (Real sigma, Real eTau, Real elc_tail_ratio, Real E0) 
  : itsSigma (sigma), itsETau (eTau), itsELCTailRatio (elc_tail_ratio), itsE0 (E0)
{
  checkInput ();
  allocateClasses ();
  return;
}

calorimeterLSF::~calorimeterLSF ()
{
  freeClasses ();
  return;
}

void calorimeterLSF::checkInput ()
{
  if (itsSigma < 0.) {
    itsSigma = 0.;
    cerr << "calorimeterLSF: sigma < 0, setting sigma = 0\n";
  }
  if (itsETau < 0.) {
    itsETau = 0.;
    cerr << "calorimeterLSF: etau < 0, setting etau = 0\n";
  }
  if (itsELCTailRatio < 0.) {
    itsELCTailRatio = 0.;
    cerr << "calorimeterLSF: nflat < 0, setting nflat = 0\n";
  }
  return;
}

void calorimeterLSF::allocateClasses ()
{
  itsELC = new electronLossContinuum (itsSigma, itsETau, itsELCTailRatio, itsE0);
  return;
}

void calorimeterLSF::freeClasses ()
{
  delete itsELC;
  return;
}

double calorimeterLSF::integrand (double deltaE) 
{
  return itsELC->getELC (deltaE);
}

Real calorimeterLSF::getLSF (Real deltaE1, Real deltaE2)
{
  //  Real answer = qags (deltaE2, deltaE1);
  //  Real answer = qng (deltaE2, deltaE1);
  setEpsRel (1.e-3); // desired accuracy of integral over energy bin
  Real limit = 1.e-10 / itsETau; // below a certain threshold throw out the tail
  if (compare (integrand (deltaE1), limit) == -1) {
    if (compare (integrand (deltaE2), limit) == -1) {
      return 0.;
    }
  }
  Real answer = qag (deltaE2, deltaE1, 4);
  int status = getStatus ();
  if (status) {
    cout << deltaE1 << " " << deltaE2 << " " << answer << "\n";
    cout << integrand (deltaE1) << " " << integrand (deltaE2) << "\n";
  }
  return qags (deltaE2, deltaE1);
}
