/***************************************************************************
    GratingProfile.cpp - ad hoc models for the line profile of a reflection 
        grating. Assumes a Gaussian core and Lorentzian wings.


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

#include <iostream>
using namespace std;

#include "xsTypes.h"
#include "Lorentzian.h"
#include "Gaussian.h"
#include "Utilities.h"
#include "isisCPPFunctionWrapper.h"

static const size_t GRATPROF_N_PARAMETERS (4);
static const size_t GRATPR2_N_PARAMETERS (6);

extern "C" void gratprof
(const RealArray& energy, const RealArray& parameter, 
   /*@unused@*/ int spectrum, RealArray& flux, 
   /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init);

extern "C" void gratpr2
(const RealArray& energy, const RealArray& parameter, 
   /*@unused@*/ int spectrum, RealArray& flux, 
   /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init);

extern "C" void C_gratprof
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init)
{
  isisCPPFunctionWrapper (energy, Nflux, parameter, spectrum, flux, 
			  fluxError, init, GRATPROF_N_PARAMETERS, &gratprof);
  return;
}

extern "C" void C_gratpr2
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init)
{
  isisCPPFunctionWrapper (energy, Nflux, parameter, spectrum, flux, 
			  fluxError, init, GRATPR2_N_PARAMETERS, &gratpr2);
  return;
}

void gratprofProcessParameter 
(const RealArray& parameter, Real& LineEnergy, Real& Sigma, Real& Gamma, 
 Real& GaussianFraction);

void gratproeProcessParameter 
(const RealArray& parameter, Real& LineEnergy, Real& Sigma, Real& Gamma, 
 Real& GaussianFraction);

void gratpr2ProcessParameter 
(const RealArray& parameter, Real& LineEnergy, Real& Sigma, Real& Gamma1, 
 Real& Gamma2, Real& GaussianFraction, Real& Ratio);

void gratpr2eProcessParameter 
(const RealArray& parameter, Real& LineEnergy, Real& Sigma, Real& Gamma1, 
 Real& Gamma2, Real& GaussianFraction, Real& Ratio);

void gratingProfile
(const RealArray& energy, RealArray& flux, Real LineEnergy, Real Sigma, 
 Real Gamma, Real GaussianFraction);

void gratingProfile2 (const RealArray& energy, RealArray& flux, 
		      Real LineEnergy, Real Sigma, 
		      Real Gamma1, Real Gamma2, 
		      Real GaussianFraction, Real Ratio);

void gratprof
(const RealArray& energy, const RealArray& parameter, 
   /*@unused@*/ int spectrum, RealArray& flux, 
   /*@unused@*/ RealArray& fluxError,
   /*@unused@*/ const string& init)
{
  fluxError.resize (0);
  flux.resize (energy.size () - 1, 0.);
  Real LineEnergy, Sigma, Gamma, GaussianFraction;
  /* Sigma and Gamma are dimensionless (relative to LineEnergy) */
  if (init == "energy") {
    gratproeProcessParameter (parameter, LineEnergy, Sigma, Gamma, 
			      GaussianFraction);
  } else {
    gratprofProcessParameter (parameter, LineEnergy, Sigma, Gamma, 
			      GaussianFraction);
  }
  gratingProfile (energy, flux, LineEnergy, Sigma, Gamma, GaussianFraction);
  return;
}

void gratpr2
(const RealArray& energy, const RealArray& parameter, 
   /*@unused@*/ int spectrum, RealArray& flux, 
 RealArray& fluxError, const string& init)
{
  fluxError.resize (0);
  flux.resize (energy.size () - 1, 0.);
  Real LineEnergy, Sigma, Gamma1, Gamma2, GaussianFraction, Ratio;
  /* Sigma and Gamma are dimensionless (relative to LineEnergy) */
  if (init == "energy") {
    gratpr2eProcessParameter (parameter, LineEnergy, Sigma, Gamma1, 
			      Gamma2, GaussianFraction, Ratio);
  } else {
    gratpr2ProcessParameter (parameter, LineEnergy, Sigma, Gamma1, 
			     Gamma2, GaussianFraction, Ratio);
  }
  gratingProfile2 (energy, flux, LineEnergy, Sigma, Gamma1, Gamma2, 
		   GaussianFraction, Ratio);
  return;
}

/***************************/

void gratingProfile
(const RealArray& energy, RealArray& flux, Real LineEnergy, Real Sigma, 
 Real Gamma, Real GaussianFraction)
{
  RealArray Lflux (flux.size (), 0.);
  RealArray Gflux (flux.size (), 0.);
  Lorentzian L (energy, LineEnergy, Gamma);
  Gaussian G (energy, LineEnergy, Sigma);
  if (G.checkOutOfBounds ()) {
    cout << "gratingProfile: Line energy (" << LineEnergy << 
      " keV) is close to or out of bounds\n";
    return;
  }
  L.getFlux (Lflux);
  G.getFlux (Gflux);
  Lflux *= (1. - GaussianFraction);
  Gflux *= GaussianFraction;
  flux = (Lflux + Gflux);
  return;
}


void gratingProfile2 (const RealArray& energy, RealArray& flux, 
		      Real LineEnergy, Real Sigma, 
		      Real Gamma1, Real Gamma2, 
		      Real GaussianFraction, Real Ratio)
{

  RealArray L1flux (flux.size (), 0.);
  RealArray L2flux (flux.size (), 0.);
  RealArray Gflux (flux.size (), 0.);

  Lorentzian L1 (energy, LineEnergy, Gamma1);
  Lorentzian L2 (energy, LineEnergy, Gamma2);
  Gaussian G (energy, LineEnergy, Sigma);
  if (G.checkOutOfBounds ()) {
    cout << "gratingProfile2: Line energy (" << LineEnergy << 
      " keV) is outside of the response range. \n";
    return;
  }
  L1.getFlux (L1flux);
  L2.getFlux (L2flux);
  G.getFlux (Gflux);
  Gflux *= GaussianFraction;
  L1flux *= Ratio * (1. - GaussianFraction) / (1. + Ratio);
  L2flux *= (1. - GaussianFraction) / (1. + Ratio);
  flux = (L1flux + L2flux + Gflux);
  return;
}


/***************************/

void gratprofProcessParameter 
(const RealArray& parameter, Real& LineEnergy, Real& Sigma, Real& Gamma, 
 Real& GaussianFraction)
{
  size_t i (0);
  Real SigmaMA = parameter[i++];
  Real GammaMA = parameter[i++];
  GaussianFraction = parameter[i++];
  Real WavelengthA = parameter[i++];
  Sigma = SigmaMA / (WavelengthA * 1.e3);
  Gamma = GammaMA / (WavelengthA * 1.e3);
  LineEnergy = convert_A_keV (WavelengthA);
  return;
}

void gratproeProcessParameter 
(const RealArray& parameter, Real& LineEnergy, Real& Sigma, Real& Gamma, 
 Real& GaussianFraction)
{
  size_t i (0);
  Real SigmaKEV = parameter[i++];
  Real GammaKEV = parameter[i++];
  GaussianFraction = parameter[i++];
  LineEnergy = parameter[i++];
  Sigma = SigmaKEV / LineEnergy;
  Gamma = GammaKEV / LineEnergy;
  return;
}

void gratpr2ProcessParameter 
(const RealArray& parameter, Real& LineEnergy, Real& Sigma, Real& Gamma1, 
 Real& Gamma2, Real& GaussianFraction, Real& Ratio)
{
  size_t i (0);
  Real SigmaMA = parameter[i++];
  Real Gamma1MA = parameter[i++];
  Real Gamma2MA = parameter[i++];
  GaussianFraction = parameter[i++];
  Ratio = parameter[i++];
  Real WavelengthA = parameter[i++];
  Sigma = SigmaMA / (WavelengthA * 1.e3);
  Gamma1 = Gamma1MA / (WavelengthA * 1.e3);
  Gamma2 = Gamma2MA / (WavelengthA * 1.e3);
  LineEnergy = convert_A_keV (WavelengthA);
  return;
}

void gratpr2eProcessParameter 
(const RealArray& parameter, Real& LineEnergy, Real& Sigma, Real& Gamma1, 
 Real& Gamma2, Real& GaussianFraction, Real& Ratio)
{
  size_t i (0);
  Real SigmaKEV = parameter[i++];
  Real Gamma1KEV = parameter[i++];
  Real Gamma2KEV = parameter[i++];
  GaussianFraction = parameter[i++];
  Ratio = parameter[i++];
  LineEnergy = parameter[i++];
  Sigma = SigmaKEV / LineEnergy;
  Gamma1 = Gamma1KEV / LineEnergy;
  Gamma2 = Gamma2KEV / LineEnergy;
  return;
}
