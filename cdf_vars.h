/*

    CHAOS: cdf_vars.h

    Copyright (C) 2022  Johnathan K Burchill

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

// Adapted from SLIDEM Processor

#ifndef CDF_VARS_H 
#define CDF_VARS_H

#include <stdint.h>

#include <cdf.h>

CDFstatus createVarFrom1DVar(CDFid id, char *name, long dataType, long startIndex, long stopIndex, void *buffer);
CDFstatus createVarFrom2DVar(CDFid id, char *name, long dataType, long startIndex, long stopIndex, void *buffer1D, uint8_t dimSize);

#endif // CDF_VARS_H
