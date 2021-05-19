/***************************************************************************
    LoadWindAbsorptionTables.h - load data for model windtabs from FITS tables.

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

#ifndef LOAD_WIND_ABSORPTION_TABLES
#define LOAD_WIND_ABSORPTION_TABLES

#include "xsTypes.h"

// Singleton
class KappaData
{
 public:
  static KappaData& instance ();
  KappaData ();
  const RealArray getKappa ();
  const RealArray getWavelength ();
  const RealArray getKappaHeII ();
  const RealArray getWavelengthHeII ();
  const RealArray getKappaVV (RealArray RelativeAbundances);
  const RealArray getEnergyVV ();
  Real getMu ();
  Real getMuHeII ();
  void refreshData ();
  bool checkStatus ();
  bool checkStatusHeII ();
  bool checkStatus2D ();
 private:
  bool isKappaOK;
  bool isKappaHeIIOK;
  bool isKappa2DOK;
  string itsFilename;
  string itsFilenameHeII;
  string itsFilename2D;

  RealArray itsKappa; // initialized as empty
  RealArray itsWavelength;
  RealArray itsKappaHeII; // initialized as empty
  RealArray itsWavelengthHeII;
  size_t itsAx1;
  size_t itsAx2;
  RealArray itsKappaZ; // initialized as empty
  // itsKappaZ will be a 1D representation of a 2D array
  // with dimensions itsAx1, itsAx2
  RealArray itsEnergyZ;
  const size_t itsNZ;
  size_t itsNEnergiesZ;
  RealArray itsAtomicNumber;
  RealArray itsAtomicMass;
  Real itsMu;
  Real itsMuHeII;

  void getFilenames ();
  void loadData ();
  void loadDataHeII ();
  void loadData2D ();
};

// Singleton
class TransmissionData
{
 public:
  static TransmissionData& instance ();
  const RealArray getTransmission ();
  const RealArray getTauStar ();
  //Real getMu ();
  // checkStatus will reload file if definitions changed
  // and will return true if everything if OK
  bool checkStatus ();
 private:
  TransmissionData ();
  bool isOK;
  string itsFilename;
  RealArray itsTransmission; // initialized as empty
  RealArray itsTauStar;
  void getFilename ();
  void loadData ();
};

// Singleton
class TransmissionData2D
{
 public:
  static TransmissionData2D& instance ();
  const RealArray getTransmission ();
  // getTransmission returns a 1D representation of a 2D array
  const RealArray getTauStar ();
  const RealArray getKappaRatio ();
  size_t getAx1 () {return itsAx1;}
  size_t getAx2 () {return itsAx2;}
  //Real getMu ();
  // checkStatus will reload file if definitions changed
  // and will return true if everything if OK
  bool checkStatus ();
 private:
  TransmissionData2D ();// private constructor for singleton
  bool isOK;
  string itsFilename;
  RealArray itsTransmission; // initialized as empty
  RealArray itsTauStar;
  RealArray itsKappaRatio;
  size_t itsAx1;
  size_t itsAx2;
  void getFilename ();
  void loadData ();
};



#endif//LOAD_WIND_ABSORPTION_TABLES
