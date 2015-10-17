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

int LoadKappa 
(RealArray& kappa, RealArray& kappaWavelength, Real& mu, bool HeII = false);

int LoadKappaZ (RealArray& kappa, RealArray& kappaEnergy, RealArray abundances);

int LoadTransmission
(RealArray& TransmissionTauStar, RealArray& Transmission);

int LoadTransmission2D
(RealArray& TransmissionTauStar, RealArray& TransmissionKappaRatio, 
 RealArray& Transmission2D, int& ax1, int& ax2);
/* Note that Transmission2D is a 1D RealArray representation of a 2D array with dimensions ax1 and ax2. */

#endif//LOAD_WIND_ABSORPTION_TABLES
