/***************************************************************************
    WindAbsorptionProfile.h - calculates the absorption profile of a 
    continuum absorbing resonance line (assuming complete destruction of
    scattered photons via Auger ionization), e.g. O IV 1s-2p.
    Under a reasonable set of assumptions, the profile is the cumulative
    distribution of an emission line profile.

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

#ifndef WIND_ABSORPTION_PROFILE_H
#define WIND_ABSORPTION_PROFILE_H

#include "WindProfile.h"

class WindAbsorptionProfile {
 public:
  WindAbsorptionProfile (const RealArray& energy, const RealArray& parameter);
  ~WindAbsorptionProfile ();
  void multiplyModelFlux (RealArray& flux);
 private:
  WindProfile* itsWindProfile;
  size_t itsEnergySize;
  void allocateClasses (const RealArray& energy, const RealArray& parameter);
  void freeClasses ();
};

#endif//WIND_ABSORPTION_PROFILE_H
