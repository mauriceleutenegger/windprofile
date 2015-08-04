/***************************************************************************
    LorentzianModels.cpp - contains prototypes and definitions for several 
      Gaussian xspec models. "W" refers to specification of rest wavelength,
      "V" refers to specification of shifts and widths in km/s, "C" refers to
      addition of a calibration shift parameter in mA.

                             -------------------
    begin				: October 2013
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
#include "Lorentzian.h"
#include "Utilities.h"
#include "isisCPPFunctionWrapper.h"

static const size_t WLORENTZ_N_PARAMETERS (3);
static const size_t VWLORENT_N_PARAMETERS (3);
static const size_t VWCLOREN_N_PARAMETERS (4);

extern "C" void vwlorent
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init);

extern "C" void C_vwlorent
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init);

extern "C" void vwcloren
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init);

extern "C" void C_vwcloren
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init);

extern "C" void wlorentz
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init);

extern "C" void C_wlorentz
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init);

void wlorentzProcessParameter
(const RealArray& parameter, Real& LineEnergy, Real& Gamma);

void vwlorentProcessParameter
(const RealArray& parameter, Real& LineEnergy, Real& Gamma);

void vwclorenProcessParameter
(const RealArray& parameter, Real& LineEnergy, Real& Gamma);

/*----------Function definitions below here-------------------*/

void vwlorent
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init)
{
  fluxError.resize (0);
  Real LineEnergy, Gamma;
  vwlorentProcessParameter (parameter, LineEnergy, Gamma);
  Lorentzian L (energy, LineEnergy, Gamma);
  L.getFlux (flux);
  return;
}

void vwcloren
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init)
{
  fluxError.resize (0);
  Real LineEnergy, Gamma;
  vwclorenProcessParameter (parameter, LineEnergy, Gamma);
  Lorentzian L (energy, LineEnergy, Gamma);
  L.getFlux (flux);
  return;
}

void wlorentz
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init)

{
  fluxError.resize (0);
  Real LineEnergy, Gamma;
  wlorentzProcessParameter (parameter, LineEnergy, Gamma);
  Lorentzian L (energy, LineEnergy, Gamma);
  L.getFlux (flux);
  return;
}

/*----------------Calls to isis C-C++ wrapper function--------------*/

void C_wlorentz
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init)
{
  isisCPPFunctionWrapper (energy, Nflux, parameter, spectrum, flux,
			  fluxError, init, WLORENTZ_N_PARAMETERS,
			  &wlorentz);
  return;
}

void C_vwlorent
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init)
{
  isisCPPFunctionWrapper (energy, Nflux, parameter, spectrum, flux,
			  fluxError, init, VWLORENT_N_PARAMETERS,
			  &vwlorent);
  return;
}

void C_vwcloren
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init)
{
  isisCPPFunctionWrapper (energy, Nflux, parameter, spectrum, flux,
			  fluxError, init, VWCLOREN_N_PARAMETERS,
			  &vwcloren);
  return;
}

/*-------------------Functions for processing input parameters-----------*/

void wlorentzProcessParameter
(const RealArray& parameter, Real& LineEnergy, Real& Gamma)
{
  size_t i (0);
  Real GammaMA = parameter[i++];
  Real DeltaWavelengthMA = parameter[i++];
  Real WavelengthA = parameter[i++];
  WavelengthA = ShiftWavelength (WavelengthA, 0., DeltaWavelengthMA);
  Gamma = GammaMA / (WavelengthA * 1.e3);
  LineEnergy = convert_A_keV (WavelengthA);
  return;
}

void vwlorentProcessParameter
(const RealArray& parameter, Real& LineEnergy, Real& Gamma)
{
  size_t i (0);
  Real GammaKMS = parameter[i++];
  Real DeltaVelocityKMS = parameter[i++];
  Real WavelengthA = parameter[i++];
  WavelengthA = ShiftWavelength 
    (WavelengthA, DeltaVelocityKMS, 0.);
  Gamma = convert_KMS_C (GammaKMS);
  LineEnergy = convert_A_keV (WavelengthA);
  return;
}

void vwclorenProcessParameter
(const RealArray& parameter, Real& LineEnergy, Real& Gamma)
{
  size_t i (0);
  Real GammaKMS = parameter[i++];
  Real DeltaVelocityKMS = parameter[i++];
  Real WavelengthA = parameter[i++];
  Real DeltaWavelengthMA = parameter[i++];
  WavelengthA = ShiftWavelength 
    (WavelengthA, DeltaVelocityKMS, DeltaWavelengthMA);
  Gamma = convert_KMS_C (GammaKMS);
  LineEnergy = convert_A_keV (WavelengthA);
  return;
}

