/***************************************************************************
    AngleAveragedTransmission.h - 
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

#ifndef MAL_ANGLE_AVERAGED_TRANSMISSION_H
#define MAL_ANGLE_AVERAGED_TRANSMISSION_H

#include <gsl/gsl_math.h>
#include "xsTypes.h"
#include "Utilities.h"
#include "OpticalDepth.h"
#include "mal_integration.h"


class AngleAveragedTransmission : public Integral {
 public:
  AngleAveragedTransmission (OpticalDepth* Tau);
  AngleAveragedTransmission (OpticalDepth* Tau, OpticalDepth* TauHeII, Real kappaRatio); // For radially dependent opacity.
  ~AngleAveragedTransmission ();
  Real getTransmission (Real u);
  double integrand (double mu);
  void setKappaRatio (Real kappaRatio) {itsKappaRatio = kappaRatio; return;}
  void printDebug ();
  void printNCalls ();
 private:
  Real itsU;
  Real itsPStar;
  OpticalDepth* itsOpticalDepth;
  OpticalDepth* itsOpticalDepthHeII;
  Real itsKappaRatio;
  bool isHeII;
};

#endif//MAL_ANGLE_AVERAGED_TRANSMISSION_H
