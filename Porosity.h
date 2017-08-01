/***************************************************************************
    Porosity.h   - Computes the porosity factor P = kappa_eff / kappa.

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

#ifndef POROSITY_H
#define POROSITY_H

#include <stdbool.h>
#include "xsTypes.h"
#include "Utilities.h"

class Porosity {
 public:
  Porosity 
    (Real TauClump0 = 0., bool Anisotropic = false,
     bool Prolate = false, bool Rosseland = false);
  void setParameters (Real TauClump0, bool Anisotropic,
		      bool Prolate, bool Rosseland);
  void setParameters (Real TauClump0);
  Real getPorosityFactor (Real u, Real mu);
  bool getPorous () const {return isPorous;}
 private:
  Real itsTauClump0;
  bool isAnisotropic;
  bool isProlate;
  bool isRosseland;
  bool isPorous;
  Real getTauClump (Real u);
  void checkInput ();
};

#endif//POROSITY_H
