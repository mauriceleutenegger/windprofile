/***************************************************************************
    IntegratedLuminosity.cpp
       integrates total luminosity over inverse radius 
       from 0 to u0
     

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

#include "IntegratedLuminosity.h"


IntegratedLuminosity::IntegratedLuminosity 
(Real q, Real u0, Real Umin, Velocity* V, AngleAveragedTransmission* T)
  : isTransparentCore (false), itsQ (q), itsU0 (u0), itsUmin (Umin), 
    itsVelocity (V), itsAngleAveragedTransmission (T)
{
  return;
}

IntegratedLuminosity::~IntegratedLuminosity ()
{
  return;
}

Real IntegratedLuminosity::getLuminosity ()
{
  Real L =  qag (itsUmin, itsU0);
  return L;
}

double IntegratedLuminosity::integrand (double u)
{
  Real w = itsVelocity->getVelocity (u);
  Real T = 1.; // allow for transparent core to find intrinsic Lx
  if (!isTransparentCore) {
    T = itsAngleAveragedTransmission->getTransmission (u);
    int status = itsAngleAveragedTransmission->getStatus ();
    if (status) {
      cout << "IntegratedLuminosity: AngleAveragedTransmission returned status code " << status << "\n";
      cout << "u = " << u << "\n";
      itsAngleAveragedTransmission->printDebug ();
    }
  }
  Real f = pow (u, itsQ);
  Real I = f * T / (w * w);
  return I;
}
