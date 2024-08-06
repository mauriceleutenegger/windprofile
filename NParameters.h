/***************************************************************************
    NParameters.h   - This header file hardcodes the number of parameters
                      supplied by the xspec models in lmodel.dat

                             -------------------
    begin				: June 2024
    copyright			: (C) 2024 by Maurice Leutenegger
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


#ifndef WP_NPARAMETERS_H
#define WP_NPARAMETERS_H

static const size_t WINDPROF_N_PARAMETERS (18);
static const size_t HWIND_N_PARAMETERS (18);
static const size_t HEWIND_N_PARAMETERS (21);
static const size_t ABSWIND_N_PARAMETERS (7);
static const size_t RADWIND_N_PARAMETERS (10);

#endif
// WP_NPARAMETERS_H
