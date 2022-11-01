#ifndef _TRACE_H
#define _TRACE_H

#include "shc.h"

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

int trace(ChaosCoefficients *coeffs, int startingDirection, double latitude, double longitude, double alt1km, double minAltkm, double maxAltkm, double *latitude2, double *longitude2, double *altitude2, long *stepsTaken);

int force(double t, const double y[], double f[], void *data);
int internalFieldNEC(double r, double theta, double phi, ChaosCoefficients *coeffs, double *bInt);


#endif // _TRACE_H
