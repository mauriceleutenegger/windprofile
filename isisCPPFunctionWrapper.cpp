/***************************************************************************
    isisCPPFunctionWrapper.cpp   -  Provides an XSPEC C wrapper for a CPP 
                       local model, so that it can be called by isis.

                             -------------------
    begin				: January 2009
    copyright			: (C) 2009 by Maurice Leutenegger
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
#include "isisCPPFunctionWrapper.h"
#include <string>
#include <iostream>

void isisCPPFunctionWrapper
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init, size_t NParameters,
 void (*PointerToXspecCPPFunction) 
 (const RealArray&, const RealArray&, int, RealArray&, RealArray&, 
  const string&) )
{
  RealArray energyArray (energy, Nflux + 1);
  RealArray parameterArray (parameter, NParameters + 1); 
  // There's one extra parameter - normalization.
  RealArray fluxArray (flux, Nflux);
  RealArray fluxErrorArray (fluxError, Nflux);
  string initString;
  if (init != NULL) { 
      initString = init;
  }
  (*PointerToXspecCPPFunction) 
    (energyArray, parameterArray, spectrum, fluxArray, fluxErrorArray, 
     initString);
  size_t fluxSize = fluxArray.size ();
  size_t fluxErrorSize = fluxErrorArray.size ();
  if (fluxSize == 0) {
    for (size_t i = 0; i < (size_t) Nflux; i++) {
      flux[i] = 0;
      fluxError[i] = 0;
    }
    return;
  }
  bool useError = true;
  if (fluxErrorSize == 0) {
    useError = false;
  }
  if (fluxSize < (size_t) Nflux) {
    std::cerr << "isisCPPFunctionWrapper: unexpected array " <<
      "size\nfluxArray.size = " << fluxSize << "\nNFlux = " <<
      Nflux << std::endl;
    for (size_t i = 0; i < (size_t) Nflux; i++) {
      flux[i] = 0;
      fluxError[i] = 0;
    }
    return;
  }
  for (size_t i = 0; i < (size_t) Nflux; i++) {
    flux[i] = fluxArray[i];
    //if (fluxErrorSize == (size_t) Nflux) {
    //  fluxError[i] = fluxErrorArray[i];
    //} else {
    //  fluxError[i] = 0;
    //}
  }
  // fill error on a second loop only if required
  if (useError) {
    for (size_t i = 0; i < (size_t) Nflux; i++) {
      fluxError[i] = fluxErrorArray[i];
    }
  }
  
  return;
}
