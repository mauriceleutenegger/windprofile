/***************************************************************************
    NeLikeGaussian.cpp -  a class to get gaussian lines for all three 
                        2p-3s lines in a Ne-like ion 

                             -------------------
    begin				: May 2021
    copyright			: (C) 2021 by Maurice Leutenegger
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

#include "NeLikeGaussian.h"
#include "Gaussian.h"
#include "Utilities.h"
#include "AtomicParameters.h"
#include <iostream>
using namespace std;

NeLikeGaussian::NeLikeGaussian 
(const RealArray& energy, const RealArray& parameter) :
  itsEnergyArray (energy), its3F (0.3), itsM2_3G (1.), itsSigma (0.), 
  itsEnergy (Nlines)
  //  , itsXFraction (0.), isXFinite (false)
{
  extractParameters (parameter);
  return;
}

NeLikeGaussian::~NeLikeGaussian ()
{
  return;
}

void NeLikeGaussian::extractParameters (const RealArray& parameter)
{
  size_t i (0);
  its3F = parameter[i++];
  itsM2_3G = parameter[i++];
  Real SigmaVelocityKMS = parameter[i++];
  Real DeltaVelocityKMS = parameter[i++];
  size_t Z = size_t (parameter[i++]);
  Real DeltaWavelengthMA = parameter[i++];
  Real DeltaWavelength3FMA = parameter[i++];
  itsSigma = convert_KMS_C (SigmaVelocityKMS);

  NeLikeParameters Ne (Z);
  RealArray LineWavelength (Nlines);
  Ne.getWavelength (LineWavelength);
  for (i = 0; i < Nlines; i++) {
    Real DWMA = DeltaWavelengthMA;
    if (i == 0) { // correct wavelength of 3F
      DWMA += DeltaWavelength3FMA;
    }
    Real Wavelength = ShiftWavelength 
      (LineWavelength[i], DeltaVelocityKMS, DWMA);
    itsEnergy[i] = convert_A_keV (Wavelength);
  }
  return;
}

void NeLikeGaussian::getFlux (RealArray& flux)
{
  size_t Nflux = itsEnergyArray.size () - 1;
  flux.resize (Nflux);
  MatrixValue neflux (Nlines); 
  /* MatrixValue is from xsTypes.h - should work like an array of RealArrays */
  Gaussian G (itsEnergyArray, itsEnergy[0], itsSigma);
  for (size_t i = 0; i < Nlines; i++) {
    neflux[i].resize (Nflux);
    G.setParameters (itsEnergy[i], itsSigma);
    if (G.checkOutOfBounds ()) {
      cout << "HeLikeGaussian:getFlux () - one or more lines have "
	   << "energies outside the response range; "
	   << "setting model flux to zero.\n";
      return;
    }
    G.getFlux (neflux[i]);
  }
  RealArray Normalization (Nlines);
  Normalization[0] = its3F;
  Normalization[1] = (1. - its3F) / (1. + itsM2_3G);
  Normalization[2] = itsM2_3G * Normalization[1];
  for (size_t i = 0; i < Nlines; i++) flux += neflux[i] * Normalization [i];
  return;
}
