/***************************************************************************
    slabtabs.cpp   - XSPEC models to compute X-ray transmission of a slab
                     of ionized material, as in a stellar wind.


                             -------------------
    begin				: Summer 2008
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

#include "xsTypes.h"
#include <gsl/gsl_const_cgsm.h>
#include "Utilities.h"
#include <iostream>
#include <fstream>
#include <string>
#include "LoadWindAbsorptionTables.h"
#include "isisCPPFunctionWrapper.h"

using namespace std;

static const size_t SLABTABS_N_PARAMETERS (1);

static const Real CONST_HC_KEV_A = GSL_CONST_CGSM_PLANCKS_CONSTANT_H * 
  GSL_CONST_CGSM_SPEED_OF_LIGHT * 1.e5 / GSL_CONST_CGSM_ELECTRON_VOLT;

static const Real CONST_NH_SIGMA_CONVERSION = 1.e22 *
  GSL_CONST_CGSM_MASS_PROTON; // factor of 10^22 to convert XSPEC NH parameter

extern "C" void slabtabs
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 const string& init);

extern "C" void slamtabs
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 const string& init);


extern "C" void slabtabsisis
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init);

extern "C" void slamtabsisis
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init);

void slabAbsorption (const RealArray& energy, const RealArray& parameter, 
		     RealArray& flux, const string& init,
		     bool MassColumnDensity);

void slabtabs
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init)
{
  fluxError.resize (0);
  slabAbsorption (energy, parameter, flux, init, false);
  return;
}

void slamtabs
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init)
{
  fluxError.resize (0);
  slabAbsorption (energy, parameter, flux, init, true);
  return;
}

void slabAbsorption (const RealArray& energy, const RealArray& parameter, 
		     RealArray& flux, const string& init, 
		     bool MassColumnDensity)
{

  size_t energySize = energy.size ();
  size_t fluxSize = energySize - 1;
  flux.resize (fluxSize);
  size_t i = 0;

  Real Sigma = 0.;
  Real Column = parameter[i++]; // units depend on calling function

  // load kappa

  RealArray kappa;
  RealArray kappaWavelength;
  Real mu;
  LoadKappa (kappa, kappaWavelength, mu);

  // determine units of column parameter
  if (MassColumnDensity) {
    Sigma = Column;
  } else {
    Sigma = Column * CONST_NH_SIGMA_CONVERSION * mu;
  }

  for (i = 0; i < fluxSize; i++) {
    Real responseWavelength = 2. * CONST_HC_KEV_A / (energy[i] + energy[i+1]); 
    size_t j = BinarySearch (kappaWavelength, responseWavelength);
    Real Tau = Sigma * kappa[j];
    flux[i] = exp (-1. * Tau);
  }
  return;
}

void slabtabsisis
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init)
{
  isisCPPFunctionWrapper (energy, Nflux, parameter, spectrum, flux, 
			  fluxError, init, SLABTABS_N_PARAMETERS, &slabtabs);
  return;
}

void slamtabsisis
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init)
{
  isisCPPFunctionWrapper (energy, Nflux, parameter, spectrum, flux, 
			  fluxError, init, SLABTABS_N_PARAMETERS, &slamtabs);
  return;
}
