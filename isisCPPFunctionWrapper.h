/***************************************************************************
    isisCPPFunctionWrapper.h   -  Provides an XSPEC C wrapper for a CPP 
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

/*NParameters is the number listed in the model definition in lmodel.dat;
  this number is one smaller than the actual number of parameters in the
  array const Real* parameter, which also includes normalization. */
void isisCPPFunctionWrapper
(const Real* energy, int Nflux, const Real* parameter, int spectrum, 
 Real* flux, Real* fluxError, const char* init, size_t NParameters,
 void (*FunctionPointer) (const RealArray&, const RealArray&, int,
			  RealArray&, RealArray&, const string&));
