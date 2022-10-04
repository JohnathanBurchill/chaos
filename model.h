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
