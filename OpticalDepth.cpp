/***************************************************************************
    OpticalDepth.cpp   - Calculates the optical depth to X-rays for a stellar
                       wind along a ray p from point z.

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

#include "OpticalDepth.h"
#include <iostream>
using namespace std;

const Real OpticalDepth::MINIMUM_VELOCITY = 0.001;

/*OpticalDepth::OpticalDepth 
(Real TauStar, Real h, Real beta, bool numerical, bool anisotropic, 
 bool rosseland, bool expansion)
  : isNumerical (numerical), itsVelocity (NULL), itsPorosity (NULL),
    itsAnalyticOpticalDepth (NULL), itsNumericalOpticalDepth (NULL),
    isHeII (false)
{
  // need to save information on h if tauclump is to be varied! Or else pass it back everytime.
  if (isNumerical) {
    allocateVelocity (beta);
    allocatePorosity (TauStar, h, anisotropic, rosseland);
    allocateNumericalOpticalDepth (TauStar);
  } else {
    allocateAnalyticOpticalDepth (TauStar, h, anisotropic, expansion);
  }
  return;
  }*/

OpticalDepth::OpticalDepth 
(Real TauStar, Real h, Real beta, bool numerical, bool anisotropic, 
 bool prolate, bool rosseland, bool expansion, bool HeII)
  : isNumerical (numerical), isHeII (HeII),
    itsVelocity (NULL), itsPorosity (NULL),
    itsAnalyticOpticalDepth (NULL), itsNumericalOpticalDepth (NULL)
{
  /* 
     need to save information on h if tauclump is to be varied! 
     Or else pass it back everytime. 
  */
  /*
  if (isHeII) {
    cout << "OpticalDepth: I am using HeII\n";
  } else {
    cout << "OpticalDepth: I am not using HeII\n";
  }
  */
  if (isNumerical) {
    allocateVelocity (beta);
    allocatePorosity (TauStar, h, anisotropic, prolate, rosseland);
    allocateNumericalOpticalDepth (TauStar);
  } else {
    allocateAnalyticOpticalDepth (TauStar, h, anisotropic, expansion);
    if (isHeII) {
      cout << "OpticalDepth: isHeII is set on, but isNumerical is set off. There will be no HeII partial ionization.\n\n";
    }
  }
  return;
}

OpticalDepth::~OpticalDepth ()
{
  freeClasses ();
  return;
}

Real OpticalDepth::getOpticalDepth (Real p, Real z)
{
  if (isNumerical) {
    return itsNumericalOpticalDepth->getOpticalDepth (p, z);
  } else {
    return itsAnalyticOpticalDepth->getOpticalDepth (p, z);
  }
}

void OpticalDepth::setParameters (Real TauStar, Real h)
{
  if (isNumerical) {
    itsPorosity->setParameters (TauStar * h);
    itsNumericalOpticalDepth->setTauStar (TauStar);
  } else {
    itsAnalyticOpticalDepth->setParameters (TauStar, h);
  }
  // This version should work correctly.
}

void OpticalDepth::allocateVelocity (Real beta)
{
  itsVelocity = new Velocity (beta, MINIMUM_VELOCITY);
  return;
}

void OpticalDepth::allocatePorosity 
(Real TauStar, Real h, bool anisotropic, bool prolate, bool rosseland)
{
  Real TauClump0 = TauStar * h;
  itsPorosity = new Porosity (TauClump0, anisotropic, prolate, rosseland);
  return;
}

void OpticalDepth::allocateNumericalOpticalDepth (Real TauStar)
{
  itsNumericalOpticalDepth = new NumericalOpticalDepth 
    (TauStar, isHeII, itsPorosity, itsVelocity);
  return;
}
void OpticalDepth::allocateAnalyticOpticalDepth 
    (Real TauStar, Real h, bool anisotropic, bool expansion)
{
  itsAnalyticOpticalDepth = new AnalyticOpticalDepth 
    (TauStar, h, anisotropic, expansion);
  return;
}

void OpticalDepth::freeClasses ()
{
  if (isNumerical) {
    delete itsNumericalOpticalDepth;
    itsNumericalOpticalDepth = NULL;
    delete itsPorosity;
    itsPorosity = NULL;
    delete itsVelocity;
    itsVelocity = NULL;
  } else {
    delete itsAnalyticOpticalDepth;
    itsAnalyticOpticalDepth = NULL;
  }
  return;
}
