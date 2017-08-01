/***************************************************************************
    WindParameter.cpp   - This class digests the input parameters passed by 
                        XSPEC and uses them to initialize the classes used
                        in the calculation of a line profile.

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

#include "WindParameter.h"
#include <iostream>
#include <gsl/gsl_const_cgsm.h>

static const size_t WINDPROF_N_PARAMETERS (18);
static const size_t HWIND_N_PARAMETERS (18);
static const size_t HEWIND_N_PARAMETERS (20);

inline Real ConvertWavelength (Real s);

// hc in kev Angstrom
const Real WindParameter::HC = 
  (1.e5 * GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_PLANCKS_CONSTANT_H 
   / GSL_CONST_CGSM_ELECTRON_VOLT); 

WindParameter::WindParameter ()
  : itsQ (0.), itsTauStar (0.), itsU0 (0.5), itsUmin (0.), 
    itsH (0.), itsTau0Star (0.), 
    itsBeta (1.), itsBetaSobolev (0.), itsKappaRatio (0.), itsR0 (1.), itsP (0.), 
    itsAtomicNumber (8),
    itsModelType (general), isNumerical (false), isAnisotropic (false),
    isProlate (false),
    isRosseland (false), isExpansion (false), isOpticallyThick (false),
    isHeII (false),
    itsWavelength (20.), itsShift (0.), itsVelocity (0.001), isVerbose (false),
    itsG (0.)
{
  return;
}

WindParameter::WindParameter (const RealArray& parameter, ModelType type)
  : itsQ (0.), itsTauStar (0.), itsU0 (0.5), itsUmin (0.), 
    itsH (0.), itsTau0Star (0.), 
    itsBeta (1.), itsBetaSobolev (0.), itsKappaRatio (0.), itsR0 (1.), itsP (0.), 
    itsAtomicNumber (8),
    itsModelType (general), isNumerical (false), isAnisotropic (false),
    isProlate (false),
    isRosseland (false), isExpansion (false), isOpticallyThick (false),
    isHeII (false),
    itsWavelength (20.), itsShift (0.), itsVelocity (0.001), isVerbose (false),
    itsG (0.)
{
  setModelType (type);
  setParameters (parameter);
  //  dump ();
  return;
}

void WindParameter::setModelType (ModelType type)
{
  itsModelType = type;
  return;
}

// Turns the array parameter into physical function parameters.
// Checks that the values are reasonable.
void WindParameter::setParameters 
(const RealArray& parameter)
{
  if (!correctNParameters (parameter.size ())) return;
  size_t i = 0;
  /*  if (itsModelType == opticaldepth) {
    setOpticalDepthParameters (parameter);
    checkInput ();
    return;
    }*/
  /*  if (itsModelType == profile) {
    setProfileParameters (parameter);
    checkInput ();
    return;
    } */
  if (itsModelType == absorption) {
    setAbsorptionParameters (parameter);
    checkInput ();
    return;
  }
  itsQ = parameter [i++];
  itsTauStar = parameter [i++];
  itsU0 = parameter [i++];
  itsUmin = parameter[i++];
  itsH = parameter [i++];
  itsTau0Star = parameter [i++];
  itsBeta = parameter [i++];
  itsBetaSobolev = parameter [i++];
  if (itsModelType == helike) {
    itsG = parameter [i++];
  }
  itsKappaRatio = parameter [i++];
  if (compare (itsKappaRatio, 0.) == 1) {
    isHeII = true;
  }
  isNumerical = bool (parameter [i++]);
  int aniso_code = int (parameter [i++]);
  if (aniso_code > 0) {
    isAnisotropic = true;
  } else {
    isAnisotropic = false;
  }
  if (aniso_code == 2) {
    isProlate = true;
  } else {
    isProlate = false;
  }
  //isAnisotropic = bool (parameter [i++]);
  isRosseland = bool (parameter [i++]);
  isExpansion = bool (parameter [i++]);
  isOpticallyThick = bool (parameter [i++]);
  if (itsModelType == general) {
    itsWavelength = parameter [i++];
  } else {
    itsAtomicNumber = int (parameter [i++]);
  }
  itsShift = ConvertWavelength (parameter [i++]);
  itsVelocity = convert_KMS_C (parameter [i++]);
  if (itsModelType == helike) {
    itsP = parameter [i++];
  }
  isVerbose = bool (parameter [i++]);
  checkInput ();
  if (itsModelType == hlike) {
    HLikeParameters HP (itsAtomicNumber);
    itsWavelength = HP.getWavelength ();
  } else if (itsModelType == helike) {
    He = HeLikeParameters (itsAtomicNumber);
    itsR0 = He.getR0 ();
    itsWavelength = 20.;
  }
  return;
}

/*
  I think that the next two functions are from the IDL
  interface, which is deprecated.
*/

/*
void WindParameter::setOpticalDepthParameters (const RealArray& parameter)
{
  size_t i (0);
  itsTauStar = parameter [i++];
  itsH = parameter [i++];
  itsBeta = parameter [i++];
  isNumerical = bool (parameter [i++]);
  isAnisotropic = bool (parameter [i++]);
  isRosseland = bool (parameter [i++]);
  isExpansion = bool (parameter [i++]);
  return;
}

void WindParameter::setProfileParameters (const RealArray& parameter)
{
  size_t i (0);
  itsQ = parameter [i++];
  itsTauStar = parameter [i++];
  itsU0 = parameter [i++];
  itsUmin = parameter[i++];
  itsH = parameter [i++];
  itsTau0Star = parameter [i++];
  itsBeta = parameter [i++];
  itsBetaSobolev = parameter [i++];
  isNumerical = bool (parameter [i++]);
  isAnisotropic = bool (parameter [i++]);
  isRosseland = bool (parameter [i++]);
  isExpansion = bool (parameter [i++]);
  isOpticallyThick = bool (parameter [i++]);
  return;
}
*/
void WindParameter::setAbsorptionParameters (const RealArray& parameter)
{
  size_t i (0);
  itsQ = parameter [i++];
  itsTauStar = parameter [i++];
  itsU0 = parameter [i++];
  itsUmin = parameter[i++];
  itsBeta = parameter [i++];
  itsWavelength = parameter [i++];
  itsVelocity = convert_KMS_C (parameter [i++]);
  return;
}


// opticaldepth and profile are passed from IDL.
// The others are passed from XSPEC and thus contain an extra parameter
// (normalization).
bool WindParameter::correctNParameters (size_t N)
{
  size_t ExpectedParameters;
  if (itsModelType == absorption) {
    ExpectedParameters = 6;
    /*  } else if (itsModelType == opticaldepth) {
    ExpectedParameters = 7;
  } else if (itsModelType == profile) {
  ExpectedParameters = 12; */
  } else if (itsModelType == helike) {
    ExpectedParameters = HEWIND_N_PARAMETERS + 1;
  } else {
    ExpectedParameters = WINDPROF_N_PARAMETERS + 1;
  }
  if (N == ExpectedParameters){
    return true;
  } else {
    cerr << "Wrong number of parameters N = " << N << ", expected " << 
      ExpectedParameters << "\n";
    return false;
  }
}

void WindParameter::checkInput ()
{
  if (compare (itsQ, -1.) != 1) {
    dump ();
    cerr << "Parameter check: Invalid q.\n";
    itsQ = 0.;
  }
  if (compare (itsTauStar, 0.) == -1) {
    dump ();
    cerr << "Parameter check: Invalid TauStar.\n";
    itsTauStar = 0.;
  }
  if ((compare (itsU0, 0.01) == -1) || (compare (itsU0, 0.99) == 1)) {
    dump ();
    cerr << "Parameter check: Invalid U0.\n";
    itsU0 = 0.5;
  }
  if (compare (itsH, 0.) == -1) {
    dump ();
    cerr << "Parameter check: Invalid h.\n";
    itsH = 0.;
  }
  if (compare (itsTau0Star, 0.) == -1) {
    dump ();
    cerr << "Parameter check: Invalid Tau0Star.\n";
    itsTau0Star = 0.;
  }
  if (compare (itsBeta, 0.) == -1) {
    dump ();
    cerr << "Parameter check: Invalid beta.\n";
    itsBeta = 0.;
  }
  if (compare (itsBetaSobolev, 0.) == -1) {
    dump ();
    cerr << "Parameter check: Invalid betaSobolev.\n";
    itsBetaSobolev = 0.;
  }
  if (compare (itsKappaRatio, 0.) == -1) {
    dump ();
    cerr << "Parameter check: Invalid kappaRatio.\n";
    itsKappaRatio = 0.;
  }
  if (compare (itsG, 0.) == -1) {
    dump ();
    cerr << "Parameter check: Invalid G.\n";
    itsG = 0.;
  }
  if (compare (itsWavelength, 1.) == -1) {
    dump ();
    cerr << "Parameter check: Invalid wavelength.\n";
    itsWavelength = 20.;
  }
  switch (itsAtomicNumber) {
  case 6: case 7: case 8: case 9: case 10: case 11: case 12: case 13:
  case 14: case 16: case 18: case 20: case 26:
    break;
  default:
    dump ();
    cerr << "Parameter check: Invalid Atomic Number.\n";
    itsAtomicNumber = 8;
    break;
  }
  if (compare (fabs(itsShift), 0.1) == 1) {
    cerr << "|Shift| " << itsShift << " (A) too large." <<
      " Maximum value is 100 mA. Setting to 0.";
    itsShift = 0.;
  }
  if (compare (itsVelocity, 0.0001) == -1) {
    cerr << "Velocity " << itsVelocity << " (c) too small." <<
      " Minumum value is about 300 km/s. Setting to about 3000 km/s.";
    itsVelocity = 0.001;
  }
  if (compare (itsP, 0.) == -1) {
    dump ();
    cerr << "Parameter check: Invalid phiRatio.\n";
    itsP = 0.;
  }
  return;
}

void WindParameter::dump ()
{
  cout << "Dumping list of parameters:\n";
  cout << "q " << itsQ << "\n";
  cout << "TauStar " << itsTauStar << "\n";
  cout << "U0 " << itsU0 << "\n";
  cout << "h " << itsH << "\n";
  cout << "Tau0Star " << itsTau0Star << "\n";
  cout << "beta " << itsBeta << "\n";
  cout << "betaSobolev " << itsBetaSobolev << "\n";
  cout << "kappaRatio " << itsKappaRatio << "\n";
  cout << "R0 " << itsR0 << "\n";
  cout << "P " << itsP << "\n";
  cout << "Z " << itsAtomicNumber << "\n";
  cout << "Wavelength (A) " << itsWavelength << "\n";
  cout << "Shift (A) " << itsShift << "\n";
  cout << "Velocity (c) " << itsVelocity << "\n";
  cout << "G " << itsG << "\n";
  cout << "P " << itsP << "\n";
  return;
}

inline Real ConvertWavelength (Real s)
{
  return (s / 1000.);
  // convert from mA to A
}

Real WindParameter::getWavelength (HeLikeType type) const
{
  return (He.getWavelength (type) + itsShift);
}

void WindParameter::setX 
(const RealArray& energy, RealArray& x, HeLikeType type)
{
  if (x.size () != energy.size ()) {
    x.resize (energy.size ());
  }
  Real wavelength;
  if (itsModelType == helike) {
    wavelength = getWavelength (type); 
  } else {
    wavelength = getWavelength ();
  }
  Real RestEnergy = HC / wavelength; // keV
  x = (RestEnergy / energy - 1.) / getVelocity ();
  return;
}

void WindParameter::initializeVelocity (Velocity*& V)
{
  V = new Velocity (itsBeta, 0.); 
  // finite minimum velocity is used only in optical depth
  return;
}

void WindParameter::initializePorosity (Porosity*& P)
{
  Real TauClump = itsTauStar * itsH;
  P = new Porosity (TauClump, isAnisotropic, isProlate, isRosseland);
  return;
}

void WindParameter::initializeOpticalDepth (OpticalDepth*& Tau, 
					    OpticalDepth*& TauHeII)
{
  /* Will allocate both regular and HeII optical depths, if (isHeII). */
  if (isHeII) {
    TauHeII = new OpticalDepth (itsTauStar, itsH, itsBeta, true, 
				isAnisotropic, isProlate,
				isRosseland, isExpansion, true);
  }
  Tau = new OpticalDepth (itsTauStar, itsH, itsBeta, isNumerical, 
			  isAnisotropic, isProlate,
			  isRosseland, isExpansion, false);
  return;
}

void WindParameter::initializeHeLikeRatio (HeLikeRatio*& He)
{
  He = new HeLikeRatio (itsR0, itsP);
  return;
}

void WindParameter::initializeResonanceScattering 
(ResonanceScattering*& RS, Velocity* V)
{
  RS = new ResonanceScattering 
    (itsTau0Star, itsBetaSobolev, isOpticallyThick, V);
  return;
}

void WindParameter::initializeLx 
(Lx*& lx, Velocity* V, HeLikeRatio* He, ResonanceScattering* RS, 
 OpticalDepth* Tau)
{
  lx = new Lx (itsQ, itsU0, itsUmin, itsBeta, wResonance, V, He, RS, Tau);
  return;
}

void WindParameter::initializeLx 
(Lx*& lx, Velocity* V, HeLikeRatio* He, ResonanceScattering* RS, 
 OpticalDepth* Tau, OpticalDepth* TauHeII)
{
  lx = new Lx (itsQ, itsU0, itsUmin, itsBeta, itsKappaRatio, 
	       wResonance, V, He, RS, Tau, TauHeII);
  return;
}


