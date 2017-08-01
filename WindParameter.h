/***************************************************************************
    WindParameter.h   - This class digests the input parameters passed by 
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

#ifndef WIND_PARAMETER_H
#define WIND_PARAMETER_H

#include <stdbool.h>
#include "Utilities.h"
#include "Porosity.h"
#include "OpticalDepth.h"
#include "HeLikeRatio.h"
#include "ResonanceScattering.h"
#include "Lx.h"
#include "AtomicParameters.h"

//enum ModelType {general, hlike, helike, opticaldepth, profile, absorption};
enum ModelType {general, hlike, helike, absorption};
/* I think that "opticaldepth" and "profile" were model types 
   which were only used in the deprecated IDL interface */

/* This class is used to turn an array of input numbers into a list
   of physical parameters. The initialize* functions then put the parameters
   into the appropriate classes. The get* functions get parameters external to
   Lx and its subsidiary classes, which are needed for a real line profile. */
class WindParameter
{
 public:
  WindParameter ();
  WindParameter (const RealArray& parameters, ModelType type);
  void setModelType (ModelType type);
  void setParameters (const RealArray& parameters);
  Real getVelocity () const {return itsVelocity;}
  Real getWavelength () const {return (itsWavelength + itsShift);}
  Real getWavelength (HeLikeType type) const;
  bool getVerbosity () const {return isVerbose;}
  Real getG () const {return itsG;}
  bool getNumerical () const {return isNumerical;}
  bool getHeII () {return isHeII;}
  void setX 
    (const RealArray& energy, RealArray& x, HeLikeType type = wResonance);
  void initializeVelocity (Velocity*& V);
  void initializePorosity (Porosity*& P);
  void initializeOpticalDepth (OpticalDepth*& Tau, OpticalDepth*& TauHeII);
  void initializeHeLikeRatio (HeLikeRatio*& He);
  void initializeResonanceScattering (ResonanceScattering*& RS, Velocity* V);
  void initializeLx (Lx*& lx, Velocity* V, HeLikeRatio* He, 
		     ResonanceScattering* RS, OpticalDepth* Tau);
  void initializeLx (Lx*& lx, Velocity* V, HeLikeRatio* He, 
		     ResonanceScattering* RS, OpticalDepth* Tau, 
		     OpticalDepth* TauHeII); // For He II recombination.
  void dump (); // Use this for debugging.
 private:
  static const Real HC;
  Real itsQ;
  Real itsTauStar;
  Real itsU0;
  Real itsUmin;
  Real itsH;
  Real itsTau0Star;
  Real itsBeta;
  Real itsBetaSobolev;
  Real itsKappaRatio;
  Real itsR0;
  Real itsP;
  int itsAtomicNumber;
  ModelType itsModelType;
  bool isNumerical;
  bool isAnisotropic;
  bool isProlate;
  bool isRosseland;
  bool isExpansion;
  bool isOpticallyThick;
  bool isHeII;
  Real itsWavelength;
  Real itsShift;
  Real itsVelocity; // in units of the speed of light
  bool isVerbose;
  HeLikeParameters He;
  Real itsG;
  bool correctNParameters (size_t N);
  void checkInput ();
  void setOpticalDepthParameters (const RealArray& parameter);
  void setProfileParameters (const RealArray& parameter);
  void setAbsorptionParameters (const RealArray& parameter);
};

#endif WIND_PARAMETER_H
