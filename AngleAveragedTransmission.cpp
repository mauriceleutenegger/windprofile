/***************************************************************************
    AngleAveragedTransmission.cpp - 
       for constant radius, integrates e^-tau over mu, 
       from -1 to the surface of occultation
     

                             -------------------
    begin				: July 2008
    copyright			: (C) 2008 by Maurice Leutenegger
    email				: Maurice.A.Leutenegger@nasa.gov
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

#include "AngleAveragedTransmission.h"

static const Real epsilon = 1.e-2;
/* 
This is the height (in stellar radii) below which the wind is occulted. 
 We do this to avoid squirrely behavior near the stellar surface where 
 the velocity law goes to zero.
*/

/* Normal version */
AngleAveragedTransmission::AngleAveragedTransmission 
(OpticalDepth* Tau) :
  itsU (0.5), itsPStar (1. + epsilon), itsOpticalDepth (Tau), itsOpticalDepthHeII (NULL), itsKappaRatio (0.), isHeII (false)
{
  return;
}

/* He recombination version */
AngleAveragedTransmission::AngleAveragedTransmission 
(OpticalDepth* Tau, OpticalDepth* TauHeII, Real kappaRatio) :
  itsU (0.5), itsPStar (1. + epsilon), itsOpticalDepth (Tau), itsOpticalDepthHeII (TauHeII), itsKappaRatio (kappaRatio), isHeII (true)
{
  return;
}

AngleAveragedTransmission::~AngleAveragedTransmission ()
{
  return;
}

// Breaking the integral into two pieces seems to save steps.
// Most of the shape seems to be near muOcculted.
// Factor of 1/2 in calculating T comes from mu integral:
// T = 0.5 int_{-1}^1 dmu integrand (mu)
Real AngleAveragedTransmission::getTransmission (Real u)
{
  itsU = u;
  Real muOcculted = -1. * sqrt (1. - itsPStar * itsPStar * itsU * itsU);
  //  Real alpha = 0.9;
  Real T = 0.5 * qag (muOcculted, 1.);
  //  Real T1 = qag (muOcculted, muOcculted * alpha);
  //  Real T2 = qag (muOcculted * alpha, 1.);

  //  printNCalls ();

  //  Real T = T1 +T2;
  return T;
}

double AngleAveragedTransmission::integrand (double mu)
{
  Real z = mu / itsU;
  Real p = sqrt (1. - mu * mu) / itsU;
  Real tau = itsOpticalDepth->getOpticalDepth (p, z);
  if (isHeII) {
    tau += itsKappaRatio * itsOpticalDepthHeII->getOpticalDepth (p, z);
  }
  return exp (-1. * tau);
}

void AngleAveragedTransmission::printDebug ()
{
  cout << "AngleAveragedTransmission: debug information -\n";
  cout << "itsU itsPStar itsKappaRatio\n";
  cout << itsU << "\t" << itsPStar << "\t" << itsKappaRatio << "\n";
}

void AngleAveragedTransmission::printNCalls ()
{
  size_t N = getNCalls ();
  resetNCalls ();
  cerr << "N " << N << "\n";
  return;
}
