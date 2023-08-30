/***************************************************************************
    sxsLSF.cpp   - XSPEC models to compute calorimeter LSFs for the SXS.
                     

                             -------------------
    begin				: November 2013
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
#include "Utilities.h"
#include "calorimeterLSF.h"
#include "Gaussian.h"
#include "isisCPPFunctionWrapper.h"
#include <gsl/gsl_poly.h>

static const size_t SXSLSF_N_PARAMETERS (5);

extern "C" void sxslsf
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init);

extern "C" void C_sxslsf
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init);

extern "C" void sxslsf2
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init);

extern "C" void C_sxslsf2
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init);


void exptail
(const RealArray& energy, RealArray& flux, Real sigmaKEV, Real e0KEV, Real etailKEV);

void elc
(const RealArray& energy, RealArray& flux, Real sigmaKEV, Real e0KEV);



void sxslsf
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init)
{
  size_t esize = energy.size ();
  size_t fsize = esize - 1;
  flux.resize (fsize);
  fluxError.resize (0);
  size_t i = 0;
  Real sigmaKEV = parameter[i++] / 1.e3; 
  Real e0KEV = parameter[i++] / 1.e3;
  Real etailKEV = parameter[i++] / 1.e3;
  Real ftail = parameter[i++];
  Real felc = parameter[i++];
  RealArray deltaE (energy.size ());
  deltaE = e0KEV - energy;
  Real elc_tail_ratio = felc / ftail;
  Real elc_tail_sum = felc + ftail;
  Real bin_size_KEV = energy[2] - energy[1];
  Real bin_fraction = bin_size_KEV / e0KEV;
  calorimeterLSF CLSF (sigmaKEV, etailKEV, elc_tail_ratio, e0KEV); 
  for (size_t j=0; j<fsize; j++) {
    flux[j] = CLSF.getLSF (deltaE[j], deltaE[j+1]);
  }
  Real total = flux.sum ();
  cout << "CLSF total:\t" << total << endl; 
  flux *= elc_tail_sum;
  Real sigmaC = sigmaKEV / e0KEV;
  Gaussian G (energy, e0KEV, sigmaC);
  RealArray gflux (fsize);
  G.getFlux (gflux);
  cout << "gflux total:\t" << gflux.sum () << endl; 
  flux += gflux * (1. - ftail);
  return;
}

void sxslsf2 
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init)
{
  // set up output arrays
  size_t esize = energy.size ();
  size_t fsize = esize - 1;
  flux.resize (fsize);
  fluxError.resize (0);
  size_t i = 0;
  // parse parameters
  Real sigmaKEV = parameter[i++] / 1.e3; 
  Real e0KEV = parameter[i++] / 1.e3;
  Real etailKEV = parameter[i++] / 1.e3;
  Real ftail = parameter[i++];
  Real felc = parameter[i++];
  // for convenience 
  RealArray deltaE (energy.size ());
  deltaE = e0KEV - energy;
  // will add components for three parts of the model
  // 1. make a normalized Gaussian for core
  Real sigmaC = sigmaKEV / e0KEV;
  Gaussian G (energy, e0KEV, sigmaC);
  RealArray gflux (fsize);
  G.getFlux (gflux);
  // 2. compute exp tail using approximation derived by Mike W
  RealArray tflux (fsize);
  exptail (energy, tflux, sigmaKEV, e0KEV, etailKEV);
  tflux /= tflux.sum ();
  // 3. approximate ELC with erf
  RealArray cflux (fsize);
  elc (energy, cflux, sigmaKEV, e0KEV);
  cflux /= cflux.sum ();
  // 4. add together
  // RENORMALIZE FIRST!
  flux += gflux * (1. - ftail - felc);
  flux += tflux * ftail;
  flux += cflux * felc;
  // normalize and then use fractions 
  return;
}

void C_sxslsf
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init)
{
  isisCPPFunctionWrapper (energy, Nflux, parameter, spectrum, flux, 
			  fluxError, init, SXSLSF_N_PARAMETERS, &sxslsf);
  return;
}

void C_sxslsf2
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init)
{
  isisCPPFunctionWrapper (energy, Nflux, parameter, spectrum, flux, 
			  fluxError, init, SXSLSF_N_PARAMETERS, &sxslsf2);
  return;
}

/*

2. should the additive Gaussian be incorporated into the calorimeterLSF class?
3. where does energy translate into deltaE?

*/

// analytic function for exponential tail
// to replace earlier quadrature formulation
void exptail
(const RealArray& energy, RealArray& flux,
 Real sigmaKEV, Real e0KEV, Real etailKEV)
{
  size_t esize = energy.size ();
  size_t fsize = esize - 1;
  RealArray x (esize);
  x = energy - e0KEV;
  RealArray cutpoint (esize);
  cutpoint = x * etailKEV / pow(sigmaKEV, 2);
  const Real p = 0.3275911;
  const Real a[] = {0.254829592, -0.284496736, 1.421413741, -1.453152027,
  1.061405429};
  RealArray denominator (esize);
  denominator = (etailKEV + p * abs (sigmaKEV * (1 + cutpoint) / sqrt (2)));
  RealArray tbar (esize);
  tbar = etailKEV / denominator;
  RealArray qbar (esize);
  for (size_t i=0; i<esize; i++) {
    qbar[i] = gsl_poly_eval(a, 5, tbar[i]);
  }
  RealArray gauss (esize);
  gauss = exp (-0.5 * pow (x/sigmaKEV, 2));
  RealArray g (esize);
  g = 0.5 * (qbar / denominator) * gauss;
  // it's easier to loop over the remaining steps
  // to avoid making a mask, and also avoid calculating h
  Real tail;
  Real h;
  Real previous;
  Real current;
  Real exp_s_t = exp (0.5 * pow (sigmaKEV/etailKEV, 2));
  previous = g[0];
  if (compare (cutpoint[0], -1) == -1) {
    tail = exp (x[0] / etailKEV);

    h = tail * exp_s_t / etailKEV;
    previous = h - g[0];
  }
  // at values where it could diverge
  for (size_t i=0; i<fsize; i++){
    if (compare (cutpoint[i+1], -1) == -1) {
          tail = exp (x[i+1] / etailKEV);
	  h = tail * exp_s_t / etailKEV;
	  current = h - g[i+1];
    } else {
      current = g[i+1];
    }
    flux[i] += (current + previous) / 2;
    previous = current;
  }
  return;
}

void elc
(const RealArray& energy, RealArray& flux, Real sigmaKEV, Real e0KEV)
{
  size_t esize = energy.size ();
  size_t fsize = esize - 1;
  RealArray dE (esize);
  dE = energy - e0KEV;
  RealArray x (esize);
  x = dE / sigmaKEV;
  Real previous = erfc(x[0]);
  RealArray answer (fsize);
  for (size_t i=0; i < fsize; i++) {
    Real current = erfc(x[i+1]);
    flux[i] = (current + previous) / 2.;
    previous = current;
    }
  return;
}
