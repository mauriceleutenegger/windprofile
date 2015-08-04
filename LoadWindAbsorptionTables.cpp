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

using namespace std;
using namespace CCfits;

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

