/***************************************************************************
    windprof.cpp   - XSPEC models to compute X-ray emission line profiles
                     from O star winds.

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

#include "xsTypes.h"
#include "WindProfile.h"
#include "WindAbsorptionProfile.h"
#include "isisCPPFunctionWrapper.h"

static const size_t WINDPROF_N_PARAMETERS (18);
static const size_t HWIND_N_PARAMETERS (18);
static const size_t HEWIND_N_PARAMETERS (20);

extern "C" void windprof
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init);

extern "C" void windprofisis
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init);

extern "C" void hwind
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init);

extern "C" void hwindisis
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init);

extern "C" void hewind
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init);

extern "C" void hewindisis
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init);

extern "C" void abswind
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init);

void windprof
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init)
{
  fluxError.resize (0);
  WindProfile W (energy, parameter, general);
  W.getModelFlux (flux);
  return;
}

void hwind
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init)
{
  fluxError.resize (0);
  WindProfile W (energy, parameter, hlike);
  W.getModelFlux (flux);
  return;
}

void hewind
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init)
{
  fluxError.resize (0);
  WindProfile W (energy, parameter, helike);
  W.getModelFlux (flux);
  return;
}

void abswind
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init)
{
  fluxError.resize (0);
  WindAbsorptionProfile W (energy, parameter);
  W.multiplyModelFlux (flux);
  return;
}

/*-------------------isis C-C++ wrapper functions----------------------*/

void windprofisis
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init)
{
  isisCPPFunctionWrapper (energy, Nflux, parameter, spectrum, flux, fluxError,
			  init, WINDPROF_N_PARAMETERS, &windprof);
  return;
}

void hwindisis
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init)
{
  isisCPPFunctionWrapper (energy, Nflux, parameter, spectrum, flux, fluxError,
			  init, HWIND_N_PARAMETERS, &hwind);
  return;
}

void hewindisis
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init)
{
  isisCPPFunctionWrapper (energy, Nflux, parameter, spectrum, flux, fluxError,
			  init, HEWIND_N_PARAMETERS, &hewind);
  return;
}
