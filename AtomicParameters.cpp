/***************************************************************************
    AtomicParameters.cpp   - Returns wavelengths for H and He-like ions, as
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

#include "AtomicParameters.h"
#include <iostream>

using namespace std;

AtomicParameters::AtomicParameters (int AtomicNumber) 
  : itsAtomicNumber (AtomicNumber)
{
  checkAtomicNumber ();
  return;
}

void AtomicParameters::setAtomicNumber (int AtomicNumber)
{
  itsAtomicNumber = AtomicNumber;
  checkAtomicNumber ();
  return;
}

void AtomicParameters::checkAtomicNumber ()
{
  switch (itsAtomicNumber) {
  case 6:
  case 7:
  case 8:
  case 10:
  case 11:
  case 12:
  case 13:
  case 14:
  case 16:
  case 18:
  case 20:
  case 26:
    break;
  default:
    cerr << "AtomicParameters: AtomicNumber " << itsAtomicNumber << 
      " not supported. Setting it to 8.\n";
    itsAtomicNumber = 8;
    break;
  }
  return;
}

/*----------------------------HLikeParameters--------------------------*/

HLikeParameters::HLikeParameters (int AtomicNumber) 
  : AtomicParameters (AtomicNumber), LA1Wavelength (0., NumberOfAtoms),
    LA2Wavelength (0., NumberOfAtoms)
{
  initialize ();
  return;
}

void HLikeParameters::initialize ()
{
  const Real LA1WavelengthArray[] = 
    {0., 0., 0., 0., 0., 33.7342, 24.7792, 18.9671, 0., 12.1321, 10.0232,
     8.41920, 7.17091, 6.18043, 0., 4.72735, 0., 3.73119, 0., 3.01848, 
     0., 0., 0., 0., 0., 1.77802};
  
  const Real LA2WavelengthArray[] = 
    {0., 0., 0., 0., 0., 33.7396, 24.7846, 18.9725, 0., 12.1375, 10.0286,
     8.42461, 7.17632, 6.18584, 0., 4.73276, 0., 3.73652, 0., 3.02390,
     0., 0., 0., 0., 0., 1.78344};

  LA1Wavelength = RealArray (LA1WavelengthArray, NumberOfAtoms);
  LA2Wavelength = RealArray (LA2WavelengthArray, NumberOfAtoms);  
  return;
}

Real HLikeParameters::getLyAlpha1 () const
{
  return LA1Wavelength[itsAtomicNumber - 1];
}

Real HLikeParameters::getLyAlpha2 () const
{
  return LA2Wavelength[itsAtomicNumber - 1];
}

Real HLikeParameters::getWavelength () const
{
  Real LyAlpha1 = getLyAlpha1 ();
  Real LyAlpha2 = getLyAlpha2 ();
  return (LyAlpha1 * 2. + LyAlpha2) / 3.;
}

/*-----------------------HeLikeParameters-------------------------*/

HeLikeParameters::HeLikeParameters (int AtomicNumber) 
  : AtomicParameters (AtomicNumber), WWavelength (0., NumberOfAtoms),
    XWavelength (0., NumberOfAtoms), YWavelength (0., NumberOfAtoms),
    ZWavelength (0., NumberOfAtoms), R0 (0., NumberOfAtoms),
    XFraction (0., NumberOfAtoms)
{
  initialize ();
  return;
}

void HeLikeParameters::initialize ()
{
  const Real WWavelengthArray[] = 
    {0., 0., 0., 0., 0., 40.2674, 28.7870, 21.6015, 0., 13.4473, 11.0029, 
     9.16875, 7.75730, 6.64795, 0., 5.03873, 0., 3.94907, 0., 3.17715, 
     0., 0., 0., 0., 0., 1.85040};
  
  const Real XWavelengthArray[] =
    {0., 0., 0., 0., 0., 40.7280, 29.0819, 21.8010, 0., 13.5503, 11.0802,
     9.22817, 7.80384, 6.68499, 0., 5.06314, 0., 3.96587, 0., 3.18910,
     0., 0., 0., 0., 0., 1.85541};
  
  const Real YWavelengthArray[] =
    {0., 0., 0., 0., 0., 40.7302, 29.0843, 21.8036, 0., 13.5531, 11.0832, 
     9.23121, 7.80696, 6.68819, 0., 5.06649, 0., 3.96936, 0., 3.19275, 
     0., 0., 0., 0., 0., 1.85952};
  
  const Real ZWavelengthArray[] =
    {0., 0., 0., 0., 0., 41.4718, 29.5346, 22.0974, 0., 13.6984, 11.1918,
     9.31362, 7.87212, 6.73949, 0., 5.10067, 0., 3.99415, 0., 3.21103,
     0., 0., 0., 0., 0., 1.86819};

  const Real R0Array[] = 
    {0., 0., 0., 0., 0., 11.0, 5.3, 3.7, 0., 3.1, 2.9, 2.7, 2.5, 2.3, 0., 2.04,
     0., 1.69, 0., 1.33, 0., 0., 0., 0., 0., 1.02};
  /* The values are from Porquet et al 2001, except S XV, which is from BDT72. 
     Al and Na are interpolated from P2001. Ar is interpolated from BDT72.
     The BDT values are significantly higher than the P2001 values.*/

  const Real XFractionArray[] =
    {0., 0., 0., 0., 0., 0.0, 0.0, 0.0, 0., 0.0, 0., 0.100, 0.148, 0.216, 0.,
     0.342, 0., 0.425, 0., 0.507, 0., 0., 0., 0., 0., 0.584};
  /* From BDT72. Ar is interpolated from S and Ca. */

  WWavelength = RealArray (WWavelengthArray, NumberOfAtoms);
  XWavelength = RealArray (XWavelengthArray, NumberOfAtoms);
  YWavelength = RealArray (YWavelengthArray, NumberOfAtoms);
  ZWavelength = RealArray (ZWavelengthArray, NumberOfAtoms);
  R0 = RealArray (R0Array, NumberOfAtoms);
  XFraction = RealArray (XFractionArray, NumberOfAtoms);

  return;
}

void HeLikeParameters::getWavelength (RealArray& wavelength) const
{
  size_t i (itsAtomicNumber - 1);
  wavelength.resize (4);
  wavelength[0] = WWavelength[i];
  wavelength[1] = XWavelength[i];
  wavelength[2] = YWavelength[i];
  wavelength[3] = ZWavelength[i];
  return;
}

Real HeLikeParameters::getWavelength (HeLikeType type) const
{
  size_t i (itsAtomicNumber - 1);
  switch (type) {
  case wResonance:
    return WWavelength[i];
    break;
  case xIntercombination:
    return XWavelength[i];
    break;
  case yIntercombination:
    return YWavelength[i];
    break;
  case zForbidden:
    return ZWavelength[i];
    break;
  default:
    cout << "HeLikeParameters::getWavelength: unrecognized type.\n";
    return 15.0;
    break;
  }
}

Real HeLikeParameters::getR0 () const
{
  size_t i (itsAtomicNumber - 1);
  return R0[i];
}

Real HeLikeParameters::getXFraction () const
{
  size_t i (itsAtomicNumber - 1);
  return XFraction[i];
}
