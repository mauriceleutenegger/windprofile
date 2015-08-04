/***************************************************************************
    GaussianModels.cpp - contains prototypes and definitions for several 
      Gaussian xspec models. "W" refers to specification of rest wavelength,
      "V" refers to specification of shifts and widths in km/s, "C" refers to
      addition of a calibration shift parameter in mA.

                             -------------------
    begin				: August 2007
    copyright			: (C) 2007 by Maurice Leutenegger
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
#include "Gaussian.h"
#include "HeLikeGaussian.h"
#include "AtomicParameters.h"
#include "Utilities.h"
#include "isisCPPFunctionWrapper.h"

static const size_t WGAUSS_N_PARAMETERS (3);
static const size_t VWGAUSS_N_PARAMETERS (3);
static const size_t VWCGAUSS_N_PARAMETERS (4);
static const size_t HGAUSS_N_PARAMETERS (3);
static const size_t HCGAUSS_N_PARAMETERS (4);
static const size_t HEGAUSS_N_PARAMETERS (5);
static const size_t HECGAUSS_N_PARAMETERS (6);

extern "C" void hgauss
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init);

extern "C" void C_hgauss
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init);

extern "C" void hcgauss
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init);

extern "C" void C_hcgauss
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init);

extern "C" void hegauss
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init);

extern "C" void C_hegauss
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init);

extern "C" void hecgauss
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init);

extern "C" void C_hecgauss
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init);

extern "C" void vwgauss
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init);

extern "C" void C_vwgauss
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init);

extern "C" void vwcgauss
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init);

extern "C" void C_vwcgauss
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init);

extern "C" void wgauss
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init);

extern "C" void C_wgauss
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init);

void wgaussProcessParameter
(const RealArray& parameter, Real& LineEnergy, Real& Sigma);

void vwgaussProcessParameter
(const RealArray& parameter, Real& LineEnergy, Real& Sigma);

void vwcgaussProcessParameter
(const RealArray& parameter, Real& LineEnergy, Real& Sigma);

void hgaussProcessParameter
(const RealArray& parameter, Real& LineEnergy, Real& Sigma);

void hcgaussProcessParameter
(const RealArray& parameter, Real& LineEnergy, Real& Sigma);

/*----------Function definitions below here-------------------*/


void hgauss
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init)
{
  fluxError.resize (0);
  Real LineEnergy, Sigma;
  hgaussProcessParameter (parameter, LineEnergy, Sigma);
  Gaussian G (energy, LineEnergy, Sigma);
  G.getFlux (flux);
  return;
}

void hcgauss
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init)
{
  fluxError.resize (0);
  Real LineEnergy, Sigma;
  hcgaussProcessParameter (parameter, LineEnergy, Sigma);
  Gaussian G (energy, LineEnergy, Sigma);
  G.getFlux (flux);
  return;
}

void hegauss
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init)
{
  fluxError.resize (0);
  size_t Nparameters (parameter.size ());
  RealArray NewParameter (Nparameters + 1);
  for (size_t i = 0; i < Nparameters - 1; i++)
    NewParameter[i] = parameter[i];
  NewParameter[Nparameters] = parameter[Nparameters - 1];
  /* Leave NewParameter[Nparameters-1] initialized to zero
     for zero calibration shift. */
  HeLikeGaussian He (energy, NewParameter);
  He.getFlux (flux);
  return;
}

void hecgauss
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init)
{
  fluxError.resize (0);
  HeLikeGaussian He (energy, parameter);
  He.getFlux (flux);
  return;
}

void vwgauss
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init)
{
  fluxError.resize (0);
  Real LineEnergy, Sigma;
  vwgaussProcessParameter (parameter, LineEnergy, Sigma);
  Gaussian G (energy, LineEnergy, Sigma);
  G.getFlux (flux);
  return;
}

void vwcgauss
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init)
{
  fluxError.resize (0);
  Real LineEnergy, Sigma;
  vwcgaussProcessParameter (parameter, LineEnergy, Sigma);
  Gaussian G (energy, LineEnergy, Sigma);
  G.getFlux (flux);
  return;
}

void wgauss
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init)

{
  fluxError.resize (0);
  Real LineEnergy, Sigma;
  wgaussProcessParameter (parameter, LineEnergy, Sigma);
  Gaussian G (energy, LineEnergy, Sigma);
  G.getFlux (flux);
  return;
}

/*----------------Calls to isis C-C++ wrapper function--------------*/

void C_wgauss
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init)
{
  isisCPPFunctionWrapper (energy, Nflux, parameter, spectrum, flux,
			  fluxError, init, WGAUSS_N_PARAMETERS,
			  &wgauss);
  return;
}

void C_vwgauss
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init)
{
  isisCPPFunctionWrapper (energy, Nflux, parameter, spectrum, flux,
			  fluxError, init, VWGAUSS_N_PARAMETERS,
			  &vwgauss);
  return;
}

void C_vwcgauss
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init)
{
  isisCPPFunctionWrapper (energy, Nflux, parameter, spectrum, flux,
			  fluxError, init, VWCGAUSS_N_PARAMETERS,
			  &vwcgauss);
  return;
}

void C_hgauss
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init)
{
  isisCPPFunctionWrapper (energy, Nflux, parameter, spectrum, flux,
  			  fluxError, init, HGAUSS_N_PARAMETERS,
  			  &hgauss);
  return;
}

void C_hcgauss
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init)
{
  isisCPPFunctionWrapper (energy, Nflux, parameter, spectrum, flux,
			  fluxError, init, HCGAUSS_N_PARAMETERS,
			  &hcgauss);
  return;
}

void C_hegauss
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init)
{
  isisCPPFunctionWrapper (energy, Nflux, parameter, spectrum, flux,
			  fluxError, init, HEGAUSS_N_PARAMETERS,
			  &hegauss);
  return;
}

void C_hecgauss
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init)
{
  isisCPPFunctionWrapper (energy, Nflux, parameter, spectrum, flux,
			  fluxError, init, HECGAUSS_N_PARAMETERS,
			  &hecgauss);
  return;
}

/*-------------------Functions for processing input parameters-----------*/

void wgaussProcessParameter
(const RealArray& parameter, Real& LineEnergy, Real& Sigma)
{
  size_t i (0);
  Real SigmaMA = parameter[i++];
  Real DeltaWavelengthMA = parameter[i++];
  Real WavelengthA = parameter[i++];
  WavelengthA = ShiftWavelength (WavelengthA, 0., DeltaWavelengthMA);
  Sigma = SigmaMA / (WavelengthA * 1.e3);
  LineEnergy = convert_A_keV (WavelengthA);
  return;
}

void vwgaussProcessParameter
(const RealArray& parameter, Real& LineEnergy, Real& Sigma)
{
  size_t i (0);
  Real SigmaKMS = parameter[i++];
  Real DeltaVelocityKMS = parameter[i++];
  Real WavelengthA = parameter[i++];
  WavelengthA = ShiftWavelength 
    (WavelengthA, DeltaVelocityKMS, 0.);
  Sigma = convert_KMS_C (SigmaKMS);
  LineEnergy = convert_A_keV (WavelengthA);
  return;
}

void vwcgaussProcessParameter
(const RealArray& parameter, Real& LineEnergy, Real& Sigma)
{
  size_t i (0);
  Real SigmaKMS = parameter[i++];
  Real DeltaVelocityKMS = parameter[i++];
  Real WavelengthA = parameter[i++];
  Real DeltaWavelengthMA = parameter[i++];
  WavelengthA = ShiftWavelength 
    (WavelengthA, DeltaVelocityKMS, DeltaWavelengthMA);
  Sigma = convert_KMS_C (SigmaKMS);
  LineEnergy = convert_A_keV (WavelengthA);
  return;
}

void hgaussProcessParameter
(const RealArray& parameter, Real& LineEnergy, Real& Sigma)
{
  size_t i (0);
  Real SigmaKMS = parameter[i++];
  Real DeltaVelocityKMS = parameter[i++];
  int Z = int (parameter[i++]);
  HLikeParameters H (Z);
  Real WavelengthA = H.getWavelength ();
  WavelengthA = ShiftWavelength 
    (WavelengthA, DeltaVelocityKMS, 0.);
  Sigma = convert_KMS_C (SigmaKMS);
  LineEnergy = convert_A_keV (WavelengthA);
  return;
}

void hcgaussProcessParameter
(const RealArray& parameter, Real& LineEnergy, Real& Sigma)
{
  size_t i (0);
  Real SigmaKMS = parameter[i++];
  Real DeltaVelocityKMS = parameter[i++];
  int Z = int (parameter[i++]);
  Real DeltaWavelengthMA = parameter[i++];
  HLikeParameters H (Z);
  Real WavelengthA = H.getWavelength ();
  WavelengthA = ShiftWavelength 
    (WavelengthA, DeltaVelocityKMS, DeltaWavelengthMA);
  Sigma = convert_KMS_C (SigmaKMS);
  LineEnergy = convert_A_keV (WavelengthA);
  return;
}

