/***************************************************************************
    ResonanceScattering.h   - Calculates the effects of resonance scattering 
                              on a stellar wind X-ray emission line profile.

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

#ifndef RESONANCE_SCATTERING_H
#define RESONANCE_SCATTERING_H

#include "xsTypes.h"
#include "Utilities.h"

/* getEscapeProbability returns the normalized escape probability p / <p>. */
class ResonanceScattering
{
 public:
  ResonanceScattering 
    (Real Tau0Star, Real BetaSobolev, bool OpticallyThick, Velocity* V);
  void setParameters 
    (Real Tau0Star, Real BetaSobolev, bool OpticallyThick, Velocity* V);
  Real getEscapeProbability (Real u, Real mu);
 private:
  Real itsTau0Star;
  Real itsBetaSobolev;
  bool isOpticallyThick;
  bool isOpticallyThin;
  Velocity* itsVelocity;
  Real getPAverage (Real Tau0, Real sigma);
  void checkInput ();
};

#endif//RESONANCE_SCATTERING_H
