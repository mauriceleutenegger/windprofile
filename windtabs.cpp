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
static const size_t VWINDTAB_N_PARAMETERS (14);

static const Real CONST_HC_KEV_A = GSL_CONST_CGSM_PLANCKS_CONSTANT_H * 
  GSL_CONST_CGSM_SPEED_OF_LIGHT * 1.e5 / GSL_CONST_CGSM_ELECTRON_VOLT;

extern "C" void windtabs
(const RealArray& energy, const RealArray& parameter, 
   /*@unused@*/ int spectrum, RealArray& flux, 
   /*@unused@*/ RealArray& fluxError,
   /*@unused@*/ const string& init);

extern "C" void C_windtabs
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
   Real* flux, Real* fluxError, const char* init);

extern "C" void vwindtab
(const RealArray& energy, const RealArray& parameter, 
   /*@unused@*/ int spectrum, RealArray& flux, 
   /*@unused@*/ RealArray& fluxError,
   /*@unused@*/ const string& init);

extern "C" void C_vwindtab
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
   Real* flux, Real* fluxError, const char* init);

extern "C" void vvwindta
(const RealArray& energy, const RealArray& parameter, 
   /*@unused@*/ int spectrum, RealArray& flux, 
   /*@unused@*/ RealArray& fluxError,
   /*@unused@*/ const string& init);

extern "C" void C_vvwindta
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
   Real* flux, Real* fluxError, const char* init);

void windtab1 (const RealArray& energy, RealArray& flux, Real RhoRstar);
void windtab2 (const RealArray& energy, RealArray& flux, Real RhoRstar);
void windtab3 (const RealArray& energy, RealArray& flux, Real RhoRstar,\
               RealArray abundances);
void writeKappaZ (const RealArray& kappa, const RealArray& kappaEnergy,\
                  const string& kappaOutFilename);

void vvwindta (const RealArray& energy, const RealArray& parameter, 
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
  size_t abundanceSize = 30;
  RealArray abundances (abundanceSize);
  for (size_t j=0; j<abundanceSize; j++) {
    abundances[j] = parameter[i++];
  }
  windtab3 (energy, flux, rhoRstar, abundances);
  return;
}

/* This model just fills out all implicit abundances and passes
   them through to vvwindta. */
void vwindtab (const RealArray& energy, const RealArray& parameter, 
   /*@unused@*/ int spectrum, RealArray& flux, 
   /*@unused@*/ RealArray& fluxError,
   /*@unused@*/ const string& init) 
{ 
  // Convert parameter to vvwindta format
  size_t newParameterSize = 31;
  RealArray newParameter (newParameterSize);
  /* This tells you which relative abundances are in the parameter
     array. Note that element 0 is set to True to pass Sigma*. */
  bool whichAbundances[] =
    {true,
     false, true, false, false, false, true, true, true, false, true,
     false, true, true, true, false, true, false, true, false, true,
     false, false, false, false, false, true, false, true, false, false};
  size_t j (0);
  for (size_t i=0; i<newParameterSize; i++) { 
    if (whichAbundances[i]) {
      newParameter[i] = parameter[j++];
    } else {
      newParameter[i] = 1.;
    }
  }
  vvwindta (energy, newParameter, spectrum, flux, fluxError, init);
  return;
}
  
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

void windtab3 (const RealArray& energy, RealArray& flux, Real rhoRstar,
               RealArray abundances)
{
  
  // load kappa - 2D table by Z; weight by abundances

  RealArray kappa;
  RealArray kappaEnergy;
  int status = LoadKappaZ (kappa, kappaEnergy, abundances);
  if (status) {return;}

  // Write out a file with kappa, given the abundances.
  string kappaOutFilename = getXspecVariable ("KAPPAZOUTFILE",\
                                              "kappaZ.txt");
  writeKappaZ (kappa, kappaEnergy, kappaOutFilename);
  
  // load optical depth and transmission
  // (this should be the same as for windtab1 -
  // only the kappa loading is different
  // for variable abundances)

  RealArray TransmissionTauStar;
  RealArray Transmission;
  status = LoadTransmission (TransmissionTauStar, Transmission);
  if (status) {
    cerr << "windtab3: Failed to load transmission file." << endl;
    return;
  }
  size_t fluxSize = flux.size ();
  for (size_t i = 0; i < fluxSize; i++) {
    Real centerEnergy = (energy[i] + energy[i+1]) / 2.;
    size_t j = BinarySearch (kappaEnergy, centerEnergy);
    Real TauStar = rhoRstar * kappa[j];
    size_t k = BinarySearch (TransmissionTauStar, TauStar);
    flux[i] = Transmission[k];
  }
  return;

}

// For now this writes to a text file
void writeKappaZ (const RealArray& kappa, const RealArray& kappaEnergy,\
                  const string& kappaOutFilename)
{
  ofstream fileHandle (kappaOutFilename.c_str());
  if (fileHandle.is_open())
    {
      size_t n = kappa.size ();
      for(size_t i = 0; i < n; i++){
        fileHandle << kappaEnergy[i] << "\t" << kappa[i] << "\n";
      }
      fileHandle.close();
    }
  else cout << "Unable to open file" << endl;
  return;
}


void C_windtabs
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init)
{
  isisCPPFunctionWrapper (energy, Nflux, parameter, spectrum, flux, 
			  fluxError, init, WINDTABS_N_PARAMETERS, &windtabs);
  return;
}

void C_vwindtab
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init)
{
  isisCPPFunctionWrapper (energy, Nflux, parameter, spectrum, flux, 
			  fluxError, init, VWINDTAB_N_PARAMETERS, &vwindtab);
  return;
}

