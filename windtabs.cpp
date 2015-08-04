/***************************************************************************
    windtabs.cpp   - XSPEC models to compute X-ray transmission averaged 
                     over the geometry of an O star wind.

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
#include "XspecUtilities.h"
#include <iostream>
#include <fstream>
#include <string>
#include "LoadWindAbsorptionTables.h"
#include "isisCPPFunctionWrapper.h"

using namespace std;

static const size_t WINDTABS_N_PARAMETERS (1);

static const Real CONST_HC_KEV_A = GSL_CONST_CGSM_PLANCKS_CONSTANT_H * 
  GSL_CONST_CGSM_SPEED_OF_LIGHT * 1.e5 / GSL_CONST_CGSM_ELECTRON_VOLT;

extern "C" void windtabs
(const RealArray& energy, const RealArray& parameter, 
   /*@unused@*/ int spectrum, RealArray& flux, 
   /*@unused@*/ RealArray& fluxError,
   /*@unused@*/ const string& init);

extern "C" void windtabsisis
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
   Real* flux, Real* fluxError, const char* init);

void windtab1 (const RealArray& energy, RealArray& flux, Real RhoRstar);
void windtab2 (const RealArray& energy, RealArray& flux, Real RhoRstar);

void windtabs
(const RealArray& energy, const RealArray& parameter, 
   /*@unused@*/ int spectrum, RealArray& flux, 
   /*@unused@*/ RealArray& fluxError,
   /*@unused@*/ const string& init) 
{ 
  // Initialize
  
  size_t energySize = energy.size ();
  size_t fluxSize = energySize - 1;
  fluxError.resize (0);
  flux.resize (fluxSize);
  size_t i = 0;
  Real rhoRstar = parameter[i++];

  if (getXspecVariable ("HEII", "0") == "1") {
    windtab2 (energy, flux, rhoRstar);
    return;
  } else {
    windtab1 (energy, flux, rhoRstar);
    return;
  }
}

/* Standard windtabs without variable He II ionziation fraction. */
void windtab1 (const RealArray& energy, RealArray& flux, Real rhoRstar)
{
  
  // load kappa

  RealArray kappa;
  RealArray kappaWavelength;
  Real mu;
  int status = LoadKappa (kappa, kappaWavelength, mu);
  if (status) {return;}

  // load optical depth and transmission

  RealArray TransmissionTauStar;
  RealArray Transmission;
  status = LoadTransmission (TransmissionTauStar, Transmission);
  if (status) {return;}
  size_t fluxSize = flux.size ();
  for (size_t i = 0; i < fluxSize; i++) {
    Real responseWavelength = 2. * CONST_HC_KEV_A / (energy[i] + energy[i+1]); 
    size_t j = BinarySearch (kappaWavelength, responseWavelength);
    Real TauStar = rhoRstar * kappa[j];
    size_t k = BinarySearch (TransmissionTauStar, TauStar);
    flux[i] = Transmission[k];
  }
  return;
}

/* windtabs for variable HeII ionization. */
void windtab2 (const RealArray& energy, RealArray& flux, Real rhoRstar)
{
  
  /* load kappa */
 
  RealArray kappa;
  RealArray kappaWavelength;
  Real mu;
  int status = LoadKappa (kappa, kappaWavelength, mu);
  if (status) {return;}

  /* load HeII kappa */

  RealArray kappaHeII;
  RealArray kappaHeIIWavelength;
  Real temp;
  status = LoadKappa (kappaHeII, kappaHeIIWavelength, temp, true);
  if (status) {return;}

  /* Calculate kappa ratio. */

  if (kappa.size () != kappaHeII.size ()) {
    cerr << "windtab2: kappa and kappaHeII arrays are different sizes.\n";
    cerr << "kappa.size () = " << kappa.size () << "\n";
    cerr << "kappaHeII.size () = " << kappaHeII.size () << "\n";
    return;
  }
  /* Assume that the two kappas are on the same wavelength grid. */
  RealArray kappaRatio (kappa.size ());
  kappaRatio = kappaHeII / kappa;

  RealArray TransmissionTauStar;
  RealArray TransmissionKappaRatio;
  /* Note that Transmission2D is a 1D Realarray representation of a 
     2D array with size (ax1, ax2). */
  RealArray Transmission2D;
  int ax1 = 0;
  int ax2 = 0;
  LoadTransmission2D (TransmissionTauStar, TransmissionKappaRatio, 
		      Transmission2D, ax1, ax2);

  /* Now that everything is loaded:
     calculate taustar and kapparatio for each wavelength;
     look up transmission 2d for those values and assign. */
  size_t fluxSize = flux.size ();
  for (size_t i = 0; i < fluxSize; i++) {
    Real responseWavelength = 2. * CONST_HC_KEV_A / (energy[i] + energy[i+1]); 
    size_t j = BinarySearch (kappaWavelength, responseWavelength);
    Real TauStar = rhoRstar * kappa[j];
    size_t k = BinarySearch (TransmissionTauStar, TauStar);
    size_t l = BinarySearch (TransmissionKappaRatio, kappaRatio[j]);
    size_t m = k * ax1 + l;
    /* I verified the array indexing, but it would be better for it to
     be built in somehow. */
    flux[i] = Transmission2D[m]; 
  }
  return;
}

void windtabsisis
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init)
{
  isisCPPFunctionWrapper (energy, Nflux, parameter, spectrum, flux, 
			  fluxError, init, WINDTABS_N_PARAMETERS, &windtabs);
  return;
}

