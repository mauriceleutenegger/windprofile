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
#include "calorimeterLSF.h"
#include "Gaussian.h"
#include "isisCPPFunctionWrapper.h"

static const size_t SXSLSF_N_PARAMETERS (5);

extern "C" void sxslsf
(const RealArray& energy, const RealArray& parameter, 
 /*@unused@*/ int spectrum, RealArray& flux, /*@unused@*/ RealArray& fluxError,
 /*@unused@*/ const string& init);

extern "C" void C_sxslsf
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init);

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
  flux *= elc_tail_sum;
  Real sigmaC = sigmaKEV / e0KEV;
  Gaussian G (energy, e0KEV, sigmaC);
  RealArray gflux (fsize);
  G.getFlux (gflux);
  flux += gflux * (1. - ftail);
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

/*

2. should the additive Gaussian be incorporated into the calorimeterLSF class?
3. where does energy translate into deltaE?

*/
