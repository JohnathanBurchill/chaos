/*

    CHAOS: trace.h

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

#ifndef _TRACE_H
#define _TRACE_H

#include "shc.h"

#define EARTH_RADIUS_KM 6371.2

enum ChaosTraceStatus
{
    CHAOS_TRACE_OK = 0,
    CHAOS_TRACE_COEFFICIENTS,
    CHAOS_TRACE_POINTER,
    CHAOS_TRACE_GSL_ERROR
};

typedef struct TracingState
{
    ChaosCoefficients *coeffs;
    double startingDirection;
    double currentDirection;
    double speed;
} TracingState;

int initializeTracer(char *coeffDir, int year, int month, int day, ChaosCoefficients *coeffs);

int trace(ChaosCoefficients *coeffs, int startingDirection, double accuracy, double latitude, double longitude, double alt1km, double minAltkm, double maxAltkm, double *latitude2, double *longitude2, double *altitude2, long *stepsTaken);

int force(double t, const double y[], double f[], void *data);
int internalFieldNEC(double r, double theta, double phi, ChaosCoefficients *coeffs, double *bInt);


#endif // _TRACE_H
