/***************************************************************************
    AtomicParameters.h   - Returns wavelengths for H and He-like ions, as
                             well as R0 and x/x+y for He-like ions.

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

#ifndef MAL_ATOMIC_PARAMETERS_H
#define MAL_ATOMIC_PARAMETERS_H

#include "xsTypes.h"
#include "Utilities.h"

class AtomicParameters
{
 public:
  AtomicParameters (int AtomicNumber = 8);
  void setAtomicNumber (int AtomicNumber);
 protected:
  int itsAtomicNumber;
 private:
  void checkAtomicNumber ();
};

class HLikeParameters : public AtomicParameters
{
 public:
  HLikeParameters (int AtomicNumber = 8);
  Real getWavelength () const; // Get weighted average of Ly Alpha.
  Real getLyAlpha1 () const;
  Real getLyAlpha2 () const;
 private:
  static const size_t NumberOfAtoms = 26;
  RealArray LA1Wavelength;  
  RealArray LA2Wavelength;
  void initialize ();
};

class HeLikeParameters : public AtomicParameters
{
 public:
  HeLikeParameters (int AtomicNumber = 8);
  void getWavelength (RealArray& wavelength) const;
  Real getWavelength (HeLikeType type) const;
  Real getR0 () const;
  Real getXFraction () const; // Returns x / (x + y).
 private:
  static const size_t NumberOfAtoms = 26;
  RealArray WWavelength;
  RealArray XWavelength;
  RealArray YWavelength;
  RealArray ZWavelength;
  RealArray R0;
  RealArray XFraction; 
  void initialize ();
};

#endif//MAL_ATOMIC_PARAMETERS_H
