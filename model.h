/*

    CHAOS: model.c

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

#ifndef _CHAOS_MODEL_H
#define _CHAOS_MODEL_H

#include "shc.h"

#include <stdint.h>

#define EARTH_RADIUS_KM 6371.2


enum CHAOS_MODEL_STATUS
{
    CHAOS_MODEL_OK = 0,
    CHAOS_MODEL_GSL = 0,
    CHAOS_MODEL_COEFFICIENTS
};

int calculateField(double r, double theta, double phi, SHCCoefficients *coeffs, double *bn, double *be, double *bc);

int calculateResiduals(ChaosCoefficients *coeffs, int interpolationSkip, uint8_t *magVariables[], size_t nInputs, double *bCore, double *bCrust, double *dbMeas);

#endif // _CHAOS_MODEL_H
