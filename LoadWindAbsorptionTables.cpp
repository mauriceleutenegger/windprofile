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
//#include "xsFortran.h"
//#include "FunctionUtility.h"
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSFunctions/functionMap.h>

// Uncomment this for debugging output
//#define LOADWINDABSTBLS_DEBUG 1

using namespace std;
using namespace CCfits;

// ----------------- class KappaData ----------------------

KappaData& KappaData::instance ()
{
  static KappaData kappaData; // calls constructor
  return kappaData;
}

KappaData::KappaData ()
  : isKappaOK (false), isKappaHeIIOK (false), isKappa2DOK (false),
    itsAx1 (0), itsAx2 (0), itsNZ (30), itsNEnergiesZ (0),
    itsMu (1.0), itsMuHeII (1.0)
{
  getFilenames ();
  loadData ();
  loadDataHeII ();
  loadData2D ();
  return;
}

const RealArray KappaData::getKappa ()
{
  return itsKappa;
}

const RealArray KappaData::getWavelength ()
{
  return itsWavelength;
}

const RealArray KappaData::getKappaHeII ()
{
  return itsKappaHeII;
}

const RealArray KappaData::getWavelengthHeII ()
{
  return itsWavelengthHeII;
}


const RealArray KappaData::getKappaVV (RealArray RelativeAbundances,
                                       bool doHeII)
{
  RealArray MassFractions (0., itsNZ);
  RealArray kappa (0., itsNEnergiesZ);

  // calculate mass fractions from abundances and atomic masses
  Real sum = 0.;
  for (size_t i=0; i<itsNZ; i++) {
    size_t Z = i + 1;
    MassFractions[i] = FunctionUtility::getAbundance (Z) * itsAtomicMass[i] *
      RelativeAbundances[i];
    sum += MassFractions[i];
  }
  if (sum > 0.) {
    MassFractions /= sum;
  } else {
    cerr << "KappaData::getKappaVV () : sum of mass fractions is <= 0" << endl;
    return kappa; // return zeros
  }
  // if the doHeII flag is set use the mass fraction to get the opacity for HeII only
  if (doHeII) {
    for (size_t j=0; j<itsNEnergiesZ; j++) {
      //kappa[j] += itsKappaZ_HeII[j] * MassFractions[1];
      kappa = itsKappaZ_HeII * MassFractions[1]; // vector * scalar = vector
      return kappa; // only return the HeII part!
    }
  }
  // sum kappas weighted by mass fractions
  // (It would be better if there was a vectorized way to do this,
  // but I don't know what it is.)
  for (size_t i=0; i<itsNZ; i++) {
    for (size_t j=0; j<itsNEnergiesZ; j++){
      size_t k = i*itsAx1 + j; // this is the 1d representation of the 2d array
      kappa[j] += itsKappaZ[k] * MassFractions[i];
    }
  }
  return kappa;
}

const RealArray KappaData::getEnergyVV ()
{
  return itsEnergyZ;
}

Real KappaData::getMu ()
{
  return itsMu;
}

Real KappaData::getMuHeII ()
{
  return itsMuHeII;
}

// call checkStatus every time windtabs functions are called

void KappaData::refreshData ()
{
  string oldFilename = itsFilename;
  string oldFilenameHeII = itsFilenameHeII;
  string oldFilename2D = itsFilename2D;
  getFilenames ();
  if (oldFilename != itsFilename) {
    loadData ();
  }
  if (oldFilenameHeII != itsFilenameHeII) {
    loadDataHeII ();
  }
  if (oldFilename2D != itsFilename2D) {
    loadData2D ();
  }
}

bool KappaData::checkStatus ()
{
  return isKappaOK;
}

bool KappaData::checkStatusHeII ()
{
  return isKappaHeIIOK;
}

bool KappaData::checkStatus2D ()
{
  return isKappa2DOK;
}

void KappaData::getFilenames ()
{
  // Get data directory from XSPEC xset variables
  string windtabsDirectory = getXspecVariable ("WINDTABSDIRECTORY", "./");
  // Get filenames from XSPEC xset variables
  itsFilename = windtabsDirectory + "/" +
    getXspecVariable ("KAPPAFILENAME", "kappa.fits");
  itsFilenameHeII = windtabsDirectory + "/" +
    getXspecVariable ("KAPPAHEIIFILENAME", "kappaHeII.fits");
  itsFilename2D = windtabsDirectory + "/" +
     getXspecVariable ("KAPPAZFILENAME", "kappa.fits");
  // perhaps would be best to completely disable default filenames
  // to prevent dumb mistakes, but leave it for now
}

void KappaData::loadData (){
  #ifdef LOADWINDABSTBLS_DEBUG
  cout << "KappaData::loadData () running" << endl;
  #endif
  // load data from "regular" data file
  string KeywordNameMu ("mu");
  try {
    int extensionNumber (1);
    unique_ptr<FITS> pInfile 
      (new FITS (itsFilename, Read, extensionNumber, false));
    /* Note that auto_ptr is deprecated in favor of unique_ptr */
    ExtHDU& table = pInfile->currentExtension ();
    size_t NumberOfRows = table.column(1).rows ();
    table.column(1).read (itsWavelength, 1, NumberOfRows);
    table.column(2).read (itsKappa, 1, NumberOfRows);
    itsMu = 1.; // default value in case keyword doesn't exist
    try {
      table.readKey (KeywordNameMu, itsMu);
    }
    catch (FitsException& issue) {
      cerr << "KappaData::loadData: CCfits / FITSio exception:" << endl;
      cerr << issue.message () << endl;
      cerr << "(file probably doesn't have mu keyword)" << endl;
      cerr << "Using mu = 1.0" << endl;
    }
    isKappaOK = true;
  }
  catch (FitsException& issue) {
    cerr << "KappaData::loadData: CCfits / FITSio exception:" << endl;
    cerr << issue.message () << endl;
    cerr << "(file probably doesn't exist)" << endl;
    isKappaOK = false;
  }
  
  return;
}

void KappaData::loadDataHeII ()
{
  #ifdef LOADWINDABSTBLS_DEBUG
  cout << "KappaData::loadDataHeII () running" << endl;
  #endif
  string KeywordNameMu ("mu");
  try {
    int extensionNumber (1);
    unique_ptr<FITS> pInfile 
      (new FITS (itsFilenameHeII, Read, extensionNumber, false));
    ExtHDU& table = pInfile->currentExtension ();
    size_t NumberOfRows = table.column(1).rows ();
    table.column(1).read (itsWavelengthHeII, 1, NumberOfRows);
    table.column(2).read (itsKappaHeII, 1, NumberOfRows);
    itsMuHeII = 1.; // default value in case keyword doesn't exist
    try {
      table.readKey (KeywordNameMu, itsMuHeII);
    }
    catch (FitsException& issue) {
      cerr << "KappaData::loadData: CCfits / FITSio exception:" << endl;
      cerr << issue.message () << endl;
      cerr << "(file probably doesn't have mu keyword)" << endl;
      cerr << "Using mu = 1.0" << endl;
    }
    isKappaHeIIOK = true;
  }
  catch (FitsException& issue) {
    cerr << "KappaData::loadData: CCfits / FITSio exception:" << endl;
    cerr << issue.message () << endl;
    cerr << "(file probably doesn't exist)" << endl;
    isKappaHeIIOK = false;
  }
  return;
}

void KappaData::loadData2D () {
  #ifdef LOADWINDABSTBLS_DEBUG
  cout << "KappaData::loadData2D () running" << endl;
  #endif
  // Set up load of 2D array kappaZ;
  // dimensions are Z and energyx
  try {
    unique_ptr<FITS> pInfile 
      (new FITS (itsFilename2D, Read, true)); // Primary HDU - Image
    PHDU& image = pInfile->pHDU ();
    image.read (itsKappaZ); // this is a 1D representation of a 2D array
    itsAx1 = image.axis (0);
    itsAx2 = image.axis (1);
    isKappa2DOK = true;
  }
  catch (FitsException& issue) {
    cerr << "KappaData::loadData2D (): CCfits / FITSio exception:" << endl;
    cerr << issue.message () << endl;
    cerr << "(failed reading KappaZ)" << endl;
    isKappa2DOK = false;
  }
  // Load energy axis for kappa table:
  try { 
    int extensionNumber (1);
    unique_ptr<FITS> pInfile 
      (new FITS (itsFilename2D, Read, extensionNumber, false));
    ExtHDU& table = pInfile->currentExtension ();
    size_t NumberOfRows = table.column(1).rows ();
    table.column(1).read (itsEnergyZ, 1, NumberOfRows);
    itsEnergyZ *= 1.e-3; // convert from eV to keV
    itsNEnergiesZ = itsEnergyZ.size ();
    isKappa2DOK = true;
  }
  catch (FitsException& issue) {
    cerr << "KappaData::loadData2D (): CCfits / FITSio exception:" << endl;
    cerr << issue.message () << endl;
    cerr << "(file probably doesn't exist)" << endl;
    isKappa2DOK = false;
  }
  // load atomic masses:
  try {
    int extensionNumber (2);
    unique_ptr<FITS> pInfile 
      (new FITS (itsFilename2D, Read, extensionNumber, false));
    ExtHDU& table = pInfile->currentExtension ();
    size_t NumberOfRows = table.column(1).rows ();
    // change column reading
    table.column(1).read (itsAtomicNumber, 1, NumberOfRows);
    table.column(2).read (itsAtomicMass, 1, NumberOfRows);
    isKappa2DOK = true;
  }
  catch (FitsException& issue) {
    cerr << "KappaData::loadData2D (): CCfits / FITSio exception:" << endl;
    cerr << issue.message () << endl;
    cerr << "(file probably doesn't exist)" << endl;
    cerr << "File was " << itsFilename2D << endl;
    isKappa2DOK = false;
  }
  // load HeII opacity:
  try {
    int extensionNumber (3);
    unique_ptr<FITS> pInfile 
      (new FITS (itsFilename2D, Read, extensionNumber, false));
    ExtHDU& table = pInfile->currentExtension ();
    size_t NumberOfRows = table.column(1).rows ();
    // change column reading
    table.column(1).read (itsKappaZ_HeII, 1, NumberOfRows);
    isKappa2DOK = true;
  }
  catch (FitsException& issue) {
    cerr << "KappaData::loadData2D (): CCfits / FITSio exception:" << endl;
    cerr << issue.message () << endl;
    cerr << "(file probably doesn't exist)" << endl;
    cerr << "File was " << itsFilename2D << endl;
    isKappa2DOK = false;
  }
}


// ----------------- class TransmissionData ---------------

TransmissionData& TransmissionData::instance ()
{
  static TransmissionData transmissionData; // calls constructor
  return transmissionData;
}


TransmissionData::TransmissionData ()
  : isOK (false)
{
  getFilename ();
  loadData ();
  return;
}


// just return copies
// might not be the most efficient, but should be safe
// and probably not even a huge efficiency hit
const RealArray TransmissionData::getTransmission ()
{
  return itsTransmission;
}

const RealArray TransmissionData::getTauStar ()
{
  return itsTauStar;
}

// call checkStatus every time windtabs functions are called
bool TransmissionData::checkStatus ()
{
  string oldFilename = itsFilename;
  getFilename ();
  if (oldFilename != itsFilename) {
    loadData ();
  }
  return isOK;
}

void TransmissionData::getFilename ()
{
  // Get data directory from XSPEC xset variables
  string windtabsDirectory = getXspecVariable ("WINDTABSDIRECTORY", "./");
  // Get filename from XSPEC xset variables  
  itsFilename = windtabsDirectory + "/" +
    getXspecVariable ("TRANSMISSIONFILENAME", "tau_transmission_HeII.fits");
}

void TransmissionData::loadData ()
{
  #ifdef LOADWINDABSTBLS_DEBUG
  cout << "Transmission::loadTransmission running" << endl;
  #endif
  try {
    int extensionNumber (1);
    unique_ptr<FITS> pInfile 
      (new FITS (itsFilename, Read, extensionNumber, false));
    ExtHDU& table = pInfile->currentExtension ();
    size_t NumberOfRows = table.column(1).rows ();
    table.column(1).read (itsTauStar, 1, NumberOfRows);
    table.column(2).read (itsTransmission, 1, NumberOfRows);
    isOK = true;
  }
  catch (FitsException& issue) {
    cerr << "Transmission::loadTransmission: " << endl;
    cerr << "CCfits / FITSio exception:" << endl;
    cerr << issue.message () << endl;
    cerr << "(file probably doesn't exist)" << endl;
    isOK = false;
  }
  return;
}


// ----------------- class TransmissionData ---------------

TransmissionData2D& TransmissionData2D::instance ()
{
  static TransmissionData2D transmissionData2D; // calls constructor
  return transmissionData2D;
}


TransmissionData2D::TransmissionData2D ()
  : isOK (false)
{
  getFilename ();
  loadData ();
  return;
}


// just return copies
// might not be the most efficient, but should be safe
// and probably not even a huge efficiency hit
const RealArray TransmissionData2D::getTransmission ()
{
  return itsTransmission;
}

const RealArray TransmissionData2D::getTauStar ()
{
  return itsTauStar;
}

const RealArray TransmissionData2D::getKappaRatio ()
{
  return itsKappaRatio;
}

// call checkStatus every time windtabs functions are called
bool TransmissionData2D::checkStatus ()
{
  string oldFilename = itsFilename;
  getFilename ();
  if (oldFilename != itsFilename) {
    loadData ();
  }
  return isOK;
}

void TransmissionData2D::getFilename ()
{
  // Get data directory from XSPEC xset variables
  string windtabsDirectory = getXspecVariable ("WINDTABSDIRECTORY", "./");
  // Get filename from XSPEC xset variables  
  itsFilename = windtabsDirectory + "/" +
    getXspecVariable ("TRANSMISSIONFILENAME2D", "tau_transmission.fits");
  return;
}

void TransmissionData2D::loadData ()
{
  #ifdef LOADWINDABSTBLS_DEBUG
  cout << "TransmissionData2D::loadData running" << endl;
  #endif
  try {
    int extensionNumber (1);
    unique_ptr<FITS> pInfile 
      (new FITS (itsFilename, Read, extensionNumber, false));
    ExtHDU& table = pInfile->currentExtension ();
    size_t NumberOfRows = table.column(1).rows ();
    table.column(1).read (itsTauStar, 1, NumberOfRows);
    isOK = true;
  }
  catch (FitsException& issue) {
    cerr << "TransmissionData2D::loadData: " << endl;
    cerr << "CCfits / FITSio exception:" << endl;
    cerr << issue.message () << endl;
    cerr << "(file probably doesn't exist)" << endl;
    isOK = false;
  }
  if (!(isOK)) {return;}
  try {
    int extensionNumber (2);
    unique_ptr<FITS> pInfile 
      (new FITS (itsFilename, Read, extensionNumber, false));
    ExtHDU& table = pInfile->currentExtension ();
    size_t NumberOfRows = table.column(1).rows ();
    table.column(1).read (itsKappaRatio, 1, NumberOfRows);
    isOK = true;
  }
  catch (FitsException& issue) {
    cerr << "TransmissionData2D::loadData: " << endl;
    cerr << "CCfits / FITSio exception:" << endl;
    cerr << issue.message () << endl;
    cerr << "(failed reading KappaRatio)" << endl;
    isOK = false;
  }
  if (!(isOK)) {return;}
  try {
    unique_ptr<FITS> pInfile 
      (new FITS (itsFilename, Read, true)); // Primary HDU - Image
    PHDU& image = pInfile->pHDU ();
    image.read (itsTransmission); // this is a 1D representation
    itsAx1 = image.axis (0);
    itsAx2 = image.axis (1);
    isOK = true;
  }
  catch (FitsException& issue) {
    cerr << "TransmissionData2D::loadData: " << endl;
    cerr << "CCfits / FITSio exception:" << endl;
    cerr << issue.message () << endl;
    cerr << "(failed reading 2D transmission data)" << endl;
    isOK = false;
  }
  return;
}
