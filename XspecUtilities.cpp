/***************************************************************************
    XspecUtilities.cpp   - Utility functions for windprofile and Gaussians
                          which depend on Xspec include files.

                             -------------------
    begin				: October 2013
    copyright			: (C) 2013 by Maurice Leutenegger
    email				: maurice.a.leutenegger@nasa.gov
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


#include "XspecUtilities.h"
#include <XSUtil/FunctionUtils/FunctionUtility.h>

using namespace std;

string getXspecVariable (string XspecVariableName, string defaultValue) 
{
  string pvalue (FunctionUtility::getModelString (XspecVariableName));
  if (pvalue.length() && pvalue != FunctionUtility::NOT_A_KEY()) 
    return pvalue;
  return defaultValue;
}
