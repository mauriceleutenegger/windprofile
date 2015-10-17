/***************************************************************************
    LoadWindAbsorptionTables.cpp - load data for model windtabs from FITS tables.

                             -------------------
    begin				: Winter 2009
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

#include <CCfits/CCfits>
#include <iostream>
//#include "Utilities.h"
#include "LoadWindAbsorptionTables.h"
#include "XspecUtilities.h"
#include "xsFortran.h"

using namespace std;
using namespace CCfits;

int fillAbundanceArray (RealArray& ExplicitRelativeAbundances,
                        RealArray RelativeAbundances);

int getMassFractions (RealArray RelativeAbundances, RealArray& MassFractions);

int LoadKappa 
(RealArray& kappa, RealArray& kappaWavelength, Real& mu, bool HeII)
{
  // Get data directory from XSPEC xset variables
  string windtabsDirectory = getXspecVariable ("WINDTABSDIRECTORY", "./");
  // Get filename from XSPEC xset variables  
  string FITSfilename = windtabsDirectory + "/";
  if (HeII) {
    FITSfilename += getXspecVariable ("KAPPAFILENAMEHEII", "kappaHeII.fits");
  } else {
    FITSfilename += getXspecVariable ("KAPPAFILENAME", "kappa.fits");
  }
  string KeywordNameMu ("mu");

  try {
    int extensionNumber (1);
    auto_ptr<FITS> pInfile 
      (new FITS (FITSfilename, Read, extensionNumber, false));
    ExtHDU& table = pInfile->currentExtension ();
    size_t NumberOfRows = table.column(1).rows ();
    table.column(1).read (kappaWavelength, 1, NumberOfRows);
    table.column(2).read (kappa, 1, NumberOfRows);
    mu = 1.; // default value in case keyword doesn't exist
    try {
      table.readKey (KeywordNameMu, mu);
    }
    catch (FitsException& issue) {
      cerr << "LoadKappa: CCfits / FITSio exception:" << endl;
      cerr << issue.message () << endl;
      cerr << "(file probably doesn't have mu keyword)" << endl;
      cerr << "Using mu = 1.0" << endl;
    }
  }
  catch (FitsException& issue) {
    cerr << "LoadKappa: CCfits / FITSio exception:" << endl;
    cerr << issue.message () << endl;
    cerr << "(file probably doesn't exist)" << endl;
    return 1;
  }
  return 0;
}


int LoadKappaZ (RealArray& kappa, RealArray& kappaEnergy, RealArray abundances)
{
  // Get data directory from XSPEC xset variables
  string windtabsDirectory = getXspecVariable ("WINDTABSDIRECTORY", "./");
  // Get filename from XSPEC xset variables  
  string FITSfilename = windtabsDirectory + "/";
  FITSfilename += getXspecVariable ("KAPPAZFILENAME", "kappa.fits");
  //  string KeywordNameMu ("mu");

  // Set up load of 2D array kappaZ;
  // dimensions are Z and energyx
  RealArray kappaZ; // this gets appropriately resized when it is loaded
  size_t ax1 (0);
  size_t ax2 (0);
  try {
    auto_ptr<FITS> pInfile 
      (new FITS (FITSfilename, Read, true)); // Primary HDU - Image
    PHDU& image = pInfile->pHDU ();
    image.read (kappaZ); // this is a 1D representation of a 2D array
    ax1 = image.axis (0);
    ax2 = image.axis (1);
  }
  catch (FitsException& issue) {
    cerr << "LoadKappaZ: CCfits / FITSio exception:" << endl;
    cerr << issue.message () << endl;
    cerr << "(failed reading KappaZ)" << endl;
    return 1;
  }

  // Load energy axis for kappa table:
  try { 
    int extensionNumber (1);
    auto_ptr<FITS> pInfile 
      (new FITS (FITSfilename, Read, extensionNumber, false));
    ExtHDU& table = pInfile->currentExtension ();
    size_t NumberOfRows = table.column(1).rows ();
    table.column(1).read (kappaEnergy, 1, NumberOfRows);
    //    table.column(2).read (kappaZ, 1, NumberOfRows); // see if this works???
  }
  catch (FitsException& issue) {
    cerr << "LoadKappa: CCfits / FITSio exception:" << endl;
    cerr << issue.message () << endl;
    cerr << "(file probably doesn't exist)" << endl;
    return 1;
  }
  kappaEnergy *= 1.e-3; // convert from eV to keV

  // get mass fractions based on xspec abund, model abund parameter,
  // and atomic masses
  RealArray massFractions;
  getMassFractions (abundances, massFractions);


  // sum kappas weighted by mass fractions
  size_t NEnergies (kappaEnergy.size ());
  size_t NZ (massFractions.size ());
  kappa.resize (NEnergies, 0.); 
  // It would be better if there was a vectorized way to do this,
  // but I don't know what it is.
  for (size_t i=0; i<NZ; i++) {
    for (size_t j=0; j<NEnergies; j++){
      size_t k = i*ax1 + j; // this is the 1d representation of the 2d array
      kappa[j] += kappaZ[k] * massFractions[i];
    }
  }
  return 0;
}


int LoadTransmission
(RealArray& TransmissionTauStar, RealArray& Transmission) 
{
  // Get data directory from XSPEC xset variables
  string windtabsDirectory = getXspecVariable ("WINDTABSDIRECTORY", "./");
  // Get filename from XSPEC xset variables  
  string FITSfilename = windtabsDirectory + "/" +
    getXspecVariable ("TRANSMISSIONFILENAME", "tau_transmission.fits");
  try {
    int extensionNumber (1);
    auto_ptr<FITS> pInfile 
      (new FITS (FITSfilename, Read, extensionNumber, false));
    ExtHDU& table = pInfile->currentExtension ();
    size_t NumberOfRows = table.column(1).rows ();
    table.column(1).read (TransmissionTauStar, 1, NumberOfRows);
    table.column(2).read (Transmission, 1, NumberOfRows);
  }
  catch (FitsException& issue) {
    cerr << "LoadTransmission: CCfits / FITSio exception:" << endl;
    cerr << issue.message () << endl;
    cerr << "(file probably doesn't exist)" << endl;
    return 1;
  }
  return 0;
}


int LoadTransmission2D
(RealArray& TransmissionTauStar, RealArray& TransmissionKappaRatio, 
 RealArray& Transmission2D, int& ax1, int& ax2)
/* Note that Transmission2D is a 1D array that needs to be indexed in 2D */
{
  // Get data directory from XSPEC xset variables
  string windtabsDirectory = getXspecVariable ("WINDTABSDIRECTORY", "./");
  // Get filename from XSPEC xset variables  
  string FITSfilename = windtabsDirectory + "/" +
    getXspecVariable ("TRANSMISSIONFILENAME2D", "tau_transmission_HeII.fits");
  try {
    int extensionNumber (1);
    auto_ptr<FITS> pInfile 
      (new FITS (FITSfilename, Read, extensionNumber, false));
    ExtHDU& table = pInfile->currentExtension ();
    size_t NumberOfRows = table.column(1).rows ();
    table.column(1).read (TransmissionTauStar, 1, NumberOfRows);
    
  }
  catch (FitsException& issue) {
    cerr << "LoadTransmission: CCfits / FITSio exception:" << endl;
    cerr << issue.message () << endl;
    cerr << "(file probably doesn't exist; failed reading TauStar)" << endl;
    return 1;
  }
  try {
    int extensionNumber (2);
    auto_ptr<FITS> pInfile 
      (new FITS (FITSfilename, Read, extensionNumber, false));
    ExtHDU& table = pInfile->currentExtension ();
    size_t NumberOfRows = table.column(1).rows ();
    table.column(1).read (TransmissionKappaRatio, 1, NumberOfRows);
    
  }
  catch (FitsException& issue) {
    cerr << "LoadTransmission: CCfits / FITSio exception:" << endl;
    cerr << issue.message () << endl;
    cerr << "(failed reading KappaRatio)" << endl;
    return 1;
  }
  try {
    auto_ptr<FITS> pInfile 
      (new FITS (FITSfilename, Read, true)); // Primary HDU - Image
    PHDU& image = pInfile->pHDU ();
    image.read (Transmission2D); // this is a 1D representation
    ax1 = image.axis (0);
    ax2 = image.axis (1);
  }
  catch (FitsException& issue) {
    cerr << "LoadTransmission: CCfits / FITSio exception:" << endl;
    cerr << issue.message () << endl;
    cerr << "(failed reading KappaRatio)" << endl;
    return 1;
  }
  return 0;
}

int getMassFractions (RealArray RelativeAbundances, RealArray& MassFractions)
{
  size_t NElements (30);
  // should try to replace this with something better
  char* ElementNames[] =
    {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
     "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
     "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"};
  RealArray AtomicNumber;
  RealArray AtomicMass;
  RealArray ExplicitRelativeAbundances (1., NElements);
  fillAbundanceArray (ExplicitRelativeAbundances, RelativeAbundances);
  // Get data directory from XSPEC xset variables
  string windtabsDirectory = getXspecVariable ("WINDTABSDIRECTORY", "./");
  // Get filename from XSPEC xset variables  
  string FITSfilename = windtabsDirectory + "/";
  FITSfilename += getXspecVariable ("KAPPAZFILENAME", "kappaZ.fits");
  try {
    int extensionNumber (2);
    auto_ptr<FITS> pInfile 
      (new FITS (FITSfilename, Read, extensionNumber, false));
    ExtHDU& table = pInfile->currentExtension ();
    size_t NumberOfRows = table.column(1).rows ();
    // change column reading
    table.column(1).read (AtomicNumber, 1, NumberOfRows);
    table.column(2).read (AtomicMass, 1, NumberOfRows);
  }
  catch (FitsException& issue) {
    cerr << "getMassFractions: CCfits / FITSio exception:" << endl;
    cerr << issue.message () << endl;
    cerr << "(file probably doesn't exist)" << endl;
    return 1;
  }
  MassFractions.resize (NElements, 0.);
  Real sum = 0.;
  for (size_t i=0; i<NElements; i++) {
    MassFractions[i] = FGABND (ElementNames[i]) * AtomicMass[i] *
      ExplicitRelativeAbundances[i];
    sum += MassFractions[i];
  }
  // renormalize
  if (sum > 0.) {
    MassFractions /= sum;
  } else {
    cerr << "getMassFractions: sum of mass fractions is <= 0" << endl;
    return 1;
  }
  return 0;
}

int fillAbundanceArray (RealArray& ExplicitRelativeAbundances,
                        RealArray RelativeAbundances)
{
  size_t NElements (30);
  // this may not work correctly???
  bool whichAbundances[] =
    {false, true, false, false, false, true, true, true, false, true,
     false, true, true, true, false, true, false, true, false, true,
     false, false, false, false, false, true, false, true, false, false};
  size_t j (0);
  for (size_t i=0; i<NElements; i++) {
    if (whichAbundances[i]) {
        ExplicitRelativeAbundances[i] = RelativeAbundances[j];
        j++;
      }
  }
  return 0;
}
