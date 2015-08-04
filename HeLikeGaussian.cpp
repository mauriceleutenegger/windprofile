/***************************************************************************
    HeLikeGaussian.cpp -  a class to get gaussian lines for all lines in
                        a He-like triplet (including both intercombo lines).

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

#include "HeLikeGaussian.h"
#include "Gaussian.h"
#include "Utilities.h"
#include "AtomicParameters.h"
#include <iostream>
using namespace std;

HeLikeGaussian::HeLikeGaussian 
(const RealArray& energy, const RealArray& parameter) :
  itsEnergyArray (energy), itsR (1.), itsG (1.), itsSigma (0.), 
  itsEnergy (Nlines), itsXFraction (0.), isXFinite (false)
{
  extractParameters (parameter);
  return;
}

HeLikeGaussian::~HeLikeGaussian ()
{
  return;
}

void HeLikeGaussian::extractParameters (const RealArray& parameter)
{
  size_t i (0);
  itsR = parameter[i++];
  itsG = parameter[i++];
  Real SigmaVelocityKMS = parameter[i++];
  Real DeltaVelocityKMS = parameter[i++];
  size_t Z = size_t (parameter[i++]);
  Real DeltaWavelengthMA = parameter[i++];
  itsSigma = convert_KMS_C (SigmaVelocityKMS);
  HeLikeParameters He (Z);
  RealArray LineWavelength (Nlines);
  He.getWavelength (LineWavelength);
  itsXFraction = He.getXFraction ();
  isXFinite = (compare (itsXFraction, 0.) == 1);
  for (i = 0; i < Nlines; i++) {
    Real Wavelength = ShiftWavelength 
      (LineWavelength[i], DeltaVelocityKMS, DeltaWavelengthMA);
    itsEnergy[i] = convert_A_keV (Wavelength);
  }
  return;
}

void HeLikeGaussian::getFlux (RealArray& flux)
{
  size_t Nflux = itsEnergyArray.size () - 1;
  flux.resize (Nflux);
  MatrixValue heflux (Nlines); 
  /* MatrixValue is from xsTypes.h - should work like an array of RealArrays */
  Gaussian G (itsEnergyArray, itsEnergy[0], itsSigma);
  for (size_t i = 0; i < Nlines; i++) {
    heflux[i].resize (Nflux);
    G.setParameters (itsEnergy[i], itsSigma);
    if (G.checkOutOfBounds ()) {
      cout << "HeLikeGaussian:getFlux () - one or more lines have "
	   << "energies outside the response range; "
	   << "setting model flux to zero.\n";
      return;
    }
    G.getFlux (heflux[i]);
  }
  RealArray Normalization (Nlines);
  Normalization[0] = 1. / (1. + itsG);
  Real iNormalization = itsG * Normalization[0] / (1. + itsR);
  Normalization[3] = itsR * iNormalization;
  Normalization[1] = itsXFraction * iNormalization;
  Normalization[2] = (1. - itsXFraction) * iNormalization;
  for (size_t i = 0; i < Nlines; i++) flux += heflux[i] * Normalization [i];
  return;
}
