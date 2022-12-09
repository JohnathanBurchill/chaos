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

#include "model.h"
#include "shc.h"
#include "util.h"

#include <stdlib.h>
#include <time.h>
#include <signal.h>
#include <stdint.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_math.h>

extern sig_atomic_t keep_running;
extern char infoHeader[50];

int calculateField(double r, double theta, double phi, SHCCoefficients *coeffs, double *bn, double *be, double *bc)
{
	double a = EARTH_RADIUS_KM;
	double aoverr = a/r;
	double magneticDerivRN = 0.0;
	double magneticDerivThetaN = 0.0;
	double magneticDerivPhiN = 0.0;
	double br = 0.0;
	double btheta = 0.0;
	double bphi = 0.0;
	size_t gRead = 0;
	size_t hRead = 0;

	int status = 0;

    double *aoverrpowers = coeffs->aoverrpowers;
    double *derivatives = coeffs->derivatives;
    double *polynomials = coeffs->polynomials;
    int minN = coeffs->minimumN;
    int maxN = coeffs->maximumN;
    double *gnm = coeffs->gNow;
    double *hnm = coeffs->hNow;
    double g = 0.0;
    double h = 0.0;
    size_t gCoeffs = 0;
    size_t hCoeffs = 0;

	aoverrpowers[0] = aoverr * aoverr * aoverr; // For potential derivatives, (a/r)^n+2, n starting at 1
	for (int n = 1; n < maxN; n++)
	{
		aoverrpowers[n] = aoverrpowers[n-1] * aoverr;
	}

    // Derivatives are d P_l^m(cos(theta)) / d theta 
    // Exclude the Condon-Shortley phase factor
	status = gsl_sf_legendre_deriv_alt_array_e(GSL_SF_LEGENDRE_SCHMIDT, maxN, cos(theta), 1, polynomials, derivatives);
	if (status)
	{
		printf("GSL error: %s\n", gsl_strerror(status));
		return CHAOS_MODEL_GSL;
	}

    size_t lInd = 0;

	for (int n = minN; n <= maxN; n++)
	{
        g = gnm[gRead++];
        lInd = gsl_sf_legendre_array_index(n, 0);
		magneticDerivRN = g * polynomials[lInd];
		magneticDerivThetaN = g * derivatives[lInd];
		magneticDerivPhiN = 0.0;
		for (int m = 1; m <= n; m++)
		{
            lInd = gsl_sf_legendre_array_index(n, m);
            g = gnm[gRead++];
            h = hnm[hRead++];
			magneticDerivRN += (g * cos((double)m*phi) + h * sin((double)m*phi) ) * polynomials[lInd];
			magneticDerivThetaN += (g * cos((double)m*phi) + h * sin((double)m*phi) ) * derivatives[lInd];
			magneticDerivPhiN += (g * (-m*sin((double)m*phi)) + h * m*cos((double)m*phi) ) * polynomials[lInd];
		}
		magneticDerivRN *= aoverrpowers[n-1] * (-((double)n+1.0));
		magneticDerivThetaN *= aoverrpowers[n-1];
		magneticDerivPhiN *= aoverrpowers[n-1];

		// Br = -di potential by di r
		br += -magneticDerivRN;
		// Btheta = -di potential by di theta / r
		btheta += -magneticDerivThetaN;
		// Bphi = -di potential by di phi / (r*sin(theta))
		bphi += -magneticDerivPhiN;
	}

    bphi /= sin(theta);

	*bn = -btheta;
	*be = bphi;
	*bc = -br;

	return CHAOS_MODEL_OK;

}

int calculateResiduals(ChaosCoefficients *coeffs, int interpolationSkip, uint8_t *magVariables[], size_t nInputs, double *bCore, double *bCrust, double *dbMeas)
{
    int status = CHAOS_MODEL_OK;

	double degrees = M_PI / 180.0;
	double a = EARTH_RADIUS_KM;
	double inputTime = 0.;
	double r = 0.;
	double theta = 0.0 * degrees;
	double phi = 0.0 * degrees;

	double deltaT = 0.0;
	double interpolationFraction = 1.0;
	size_t lastIndex = 0;

    fprintf(stdout, "%sCalculating fields...\n", infoHeader);

    // Get fields for first time
    inputTime = ((double*)magVariables[0])[0];
    theta = (90.0 - ((double*)magVariables[1])[0]) * degrees;
    phi = ((double*)magVariables[2])[0] * degrees;
    r = ((double*)magVariables[3])[0]/1000.;

    status = calculateField(r, theta, phi, &coeffs->core, bCore, bCore+1, bCore+2);
    if (status != CHAOS_MODEL_OK)
        return status;

    status = calculateField(r, theta, phi, &coeffs->crust, bCrust, bCrust+1, bCrust+2);
    if (status != CHAOS_MODEL_OK)
        return status;

    dbMeas[0] = ((double*)magVariables[4])[0] - bCore[0] - bCrust[0];
    dbMeas[1] = ((double*)magVariables[4])[1] - bCore[1] - bCrust[1];
    dbMeas[2] = ((double*)magVariables[4])[2] - bCore[2] - bCrust[2];

    lastIndex = 0;
    for (size_t t = interpolationSkip; t < nInputs && keep_running == 1; t += interpolationSkip)
    {
        inputTime = ((double*)magVariables[0])[t];
        theta = (90.0 - ((double*)magVariables[1])[t]) * degrees;
        phi = ((double*)magVariables[2])[t] * degrees;
        r = ((double*)magVariables[3])[t]/1000.;

        status = calculateField(r, theta, phi, &coeffs->core, bCore+t*3, bCore+t*3+1, bCore+t*3+2);
        if (status != CHAOS_MODEL_OK)
            return status;

        status = calculateField(r, theta, phi, &coeffs->crust, bCrust+t*3, bCrust+t*3+1, bCrust+t*3+2);
        if (status != CHAOS_MODEL_OK)
            return status;

        dbMeas[t*3]   = ((double*)magVariables[4])[t*3 + 0] - bCore[t*3 + 0] - bCrust[t*3 + 0];
        dbMeas[t*3+1] = ((double*)magVariables[4])[t*3 + 1] - bCore[t*3 + 1] - bCrust[t*3 + 1];
        dbMeas[t*3+2] = ((double*)magVariables[4])[t*3 + 2] - bCore[t*3 + 2] - bCrust[t*3 + 2];

        // Linearly interpolate at skipped epochs
        deltaT = ((double*)magVariables[0])[t] - ((double*)magVariables[0])[t - interpolationSkip];
        if (deltaT <= 0.)
            deltaT = 1.0; // arbitrary
        for (ssize_t i = t-interpolationSkip + 1; i < t; i++)
        {
            interpolationFraction = (((double*)magVariables[0])[i] - ((double*)magVariables[0])[t-interpolationSkip]) / deltaT;
            bCore[i*3] = bCore[(t-interpolationSkip)*3] + interpolationFraction * (bCore[t*3] - bCore[(t-interpolationSkip)*3]);
            bCore[i*3+1] = bCore[(t-interpolationSkip)*3+1] + interpolationFraction * (bCore[t*3+1] - bCore[(t-interpolationSkip)*3+1]);
            bCore[i*3+2] = bCore[(t-interpolationSkip)*3+2] + interpolationFraction * (bCore[t*3+2] - bCore[(t-interpolationSkip)*3+2]);

            bCrust[i*3] = bCrust[(t-interpolationSkip)*3] + interpolationFraction * (bCrust[t*3] - bCrust[(t-interpolationSkip)*3]);
            bCrust[i*3+1] = bCrust[(t-interpolationSkip)*3+1] + interpolationFraction * (bCrust[t*3+1] - bCrust[(t-interpolationSkip)*3+1]);
            bCrust[i*3+2] = bCrust[(t-interpolationSkip)*3+2] + interpolationFraction * (bCrust[t*3+2] - bCrust[(t-interpolationSkip)*3+2]);

            // Copy skipped measured fields
            dbMeas[i*3]   = ((double*)magVariables[4])[i*3 + 0] - bCore[i*3 + 0] - bCrust[i*3 + 0];
            dbMeas[i*3+1] = ((double*)magVariables[4])[i*3 + 1] - bCore[i*3 + 1] - bCrust[i*3 + 1];
            dbMeas[i*3+2] = ((double*)magVariables[4])[i*3 + 2] - bCore[i*3 + 2] - bCrust[i*3 + 2];

        }
        lastIndex = t;

    }
    // Calculate actual field for remaining inputs (up to INTEGRATION_SKIP-1 of them)
    for (size_t t = lastIndex+1; t < nInputs && keep_running == 1; t++)
    {
        inputTime = ((double*)magVariables[0])[t];
        theta = (90.0 - ((double*)magVariables[1])[t]) * degrees;
        phi = ((double*)magVariables[2])[t] * degrees;
        r = ((double*)magVariables[3])[t]/1000.;

        status = calculateField(r, theta, phi, &coeffs->core, bCore+t*3, bCore+t*3+1, bCore+t*3+2);
        if (status != CHAOS_MODEL_OK)
            return status;

        status = calculateField(r, theta, phi, &coeffs->crust, bCrust+t*3, bCrust+t*3+1, bCrust+t*3+2);
        if (status != CHAOS_MODEL_OK)
             return status;

        dbMeas[t*3]   = ((double*)magVariables[4])[t*3 + 0] - bCore[t*3 + 0] - bCrust[t*3 + 0];
        dbMeas[t*3+1] = ((double*)magVariables[4])[t*3 + 1] - bCore[t*3 + 1] - bCrust[t*3 + 1];
        dbMeas[t*3+2] = ((double*)magVariables[4])[t*3 + 2] - bCore[t*3 + 2] - bCrust[t*3 + 2];

    }

    return CHAOS_MODEL_OK;

}
