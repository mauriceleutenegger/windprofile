#include "xsTypes.h"
#include "IntegratedLuminosity.h"
#include "AngleAveragedTransmission.h"
#include "OpticalDepth.h"
#include "Utilities.h"
#include "isisCPPFunctionWrapper.h"
#include <gsl/gsl_const_cgsm.h>
#include "LoadWindAbsorptionTables.h"

static const Real CONST_HC_KEV_A = GSL_CONST_CGSM_PLANCKS_CONSTANT_H * 
    GSL_CONST_CGSM_SPEED_OF_LIGHT * 1.e5 / GSL_CONST_CGSM_ELECTRON_VOLT;

static const size_t WINDCABS_N_PARAMETERS (5);

extern "C" void windcabs
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init);

extern "C" void C_windcabs
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init);

void windcabs
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init)
{
  // ------------------- Initialize -----------------------

  size_t energySize = energy.size ();
  size_t fluxSize = energySize - 1;
  fluxError.resize (0);
  flux.resize (fluxSize);

  // ------------------- Extract parameters -------------------

  size_t i = 0;
  Real q = parameter[i++];
  Real rhoRstar = parameter[i++];
  Real u0 = parameter[i++];
  Real umin = parameter[i++];
  Real beta = parameter[i++];

  // ------------------- Allocate classes --------------

  Velocity* V;
  OpticalDepth* Tau;
  AngleAveragedTransmission* T;
  IntegratedLuminosity* L;
  V = new Velocity (beta, 0.);
  Tau = new OpticalDepth (0., 0., beta);
  T = new AngleAveragedTransmission (Tau);
  L = new IntegratedLuminosity (q, u0, umin, V, T);
  T->setEpsRel (1.e-2);
  L->setEpsRel (1.e-2);
  Real IntrinsicLuminosity = L->getLuminosity ();

  // --------------- Load kappas -----------------

  // Singleton container class only loads when filename changes in xset
  KappaData& theKappaData = KappaData::instance ();
  // check if we need to load a new file:
  theKappaData.refreshData ();
  // make sure the data are valid
  bool Kstatus = theKappaData.checkStatus ();
  if (!Kstatus) {
    cerr << "windtab1: Problem with kappa file." << endl;
    return;
  }
  // get data from container
  RealArray kappa = theKappaData.getKappa ();
  RealArray kappaWavelength = theKappaData.getWavelength ();
  //Real mu = theKappaData.getMu ();
  // apparently I'm not using mu anywhere? so comment it out
  
  // -------------- Calculate transmission -------------------

  for (i = 0; i < fluxSize; i++) {
    Real wavelength = 2. * CONST_HC_KEV_A / (energy[i] + energy[i+1]); 
    // central wavelength
    size_t j = BinarySearch (kappaWavelength, wavelength);
    Real TauStar = rhoRstar * kappa[j];
    Tau->setParameters (TauStar, 0.);
    Real luminosity = L->getLuminosity () / IntrinsicLuminosity;
    flux[i] = luminosity;
  }

  // ---------------- Clean up -------------------------

  delete L;
  delete T;
  delete Tau;
  delete V;
  return;

}

void C_windcabs
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init)
{
  isisCPPFunctionWrapper (energy, Nflux, parameter, spectrum, flux, 
			  fluxError, init, WINDCABS_N_PARAMETERS, &windcabs);
  return;
}
