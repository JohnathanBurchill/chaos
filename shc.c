/*

    CHAOS: shc.c

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

#include "shc.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fts.h>
#include <time.h>

#include <gsl/gsl_sf_legendre.h>

int loadModelCoefficients(const char *coeffDir, ChaosCoefficients *coeffs)
{
    bzero(coeffs, sizeof(ChaosCoefficients));

	char *searchPath[2] = {NULL, NULL};
    searchPath[0] = (char *)coeffDir;

    FTS *fts = fts_open(searchPath, FTS_LOGICAL | FTS_NOSTAT, NULL);
    if (fts == NULL)
        return SHC_DIRECTORY_READ;

    FTSENT *f = fts_read(fts);

    int status = 0;

    SHCCoefficients *c = NULL;

    // This assumes only one version of CHAOS SHC files are present in the directory
    while (f)
    {
        if (f->fts_namelen >= 18 && strcmp(f->fts_name + f->fts_namelen - 8, "core.shc") == 0)
            c = &coeffs->core;
        else if (f->fts_namelen >= 31 && strcmp(f->fts_name + f->fts_namelen - 16, "extrapolated.shc") == 0)
            c = &coeffs->coreExtrapolation;
        else if (f->fts_namelen >= 20 && strcmp(f->fts_name + f->fts_namelen - 10, "static.shc") == 0)
            c = &coeffs->crust;
        else
            c = NULL;

        if (c != NULL)
        {
            strncpy((char *)c->coeffFilename, f->fts_path, FILENAME_MAX-1);
            status = loadSHCCoefficients(c);
            if (status != SHC_OK)
            {
                fts_close(fts);
                return status;
            }
        }
        f = fts_read(fts);
    }

    coeffs->initialized = coeffs->core.initialized && coeffs->coreExtrapolation.initialized && coeffs->crust.initialized;
    if (!coeffs->initialized)
        return SHC_OK;

    // Time interpolation needs each polynomial piece sampled at bSplineOrder = bSplineSteps + 1 times
    SHCCoefficients *timeDependent[2] = {&coeffs->core, &coeffs->coreExtrapolation};
    for (int i = 0; i < 2; i++)
    {
        c = timeDependent[i];
        if (c->bSplineSteps < 1 || c->bSplineOrder != c->bSplineSteps + 1 || c->bSplineOrder > SHC_MAX_SPLINE_ORDER || (c->numberOfTimes - 1) % c->bSplineSteps != 0)
        {
            fprintf(stderr, "Unexpected SHC spline order %d and step %d for %d times in %s\n", c->bSplineOrder, c->bSplineSteps, c->numberOfTimes, c->coeffFilename);
            return SHC_INTERPOLATION;
        }
    }

	// Crustal field is static, so copy g and h into gNow and hNow once
	if (coeffs->crust.numberOfTimes != 1)
	{
		fprintf(stderr, "Expected 1 static (i.e. crustal) SHC time, got %d times\n", coeffs->crust.numberOfTimes);
		return SHC_INTERPOLATION;
	}
	memcpy(coeffs->crust.gNow, coeffs->crust.gTimeSeries, coeffs->crust.gCoeffs * sizeof(double));
	memcpy(coeffs->crust.hNow, coeffs->crust.hTimeSeries, coeffs->crust.hCoeffs * sizeof(double));

    return SHC_OK;
}

int loadSHCCoefficients(SHCCoefficients *coeffs)
{
	int status = 0;
	char line[255];

	FILE *f = fopen(coeffs->coeffFilename, "r");	
	if (f == NULL)
		return SHC_FILE_READ;

	int minN = 0, maxN = 0, nTimes = 0, splineOrder = 0, splineSteps = 0; 
	size_t nTerms = 0;
    size_t nCoeffs = 0;
    size_t nInfoChars = 0;
    size_t len = 0;
	while (fgets(line, 255, f) != NULL && line[0] == '#')
    {
        len = strlen(line);
        if (nInfoChars + len < SHC_INFO_BUFFER_SIZE - 1)
        {
            strcpy((char *)coeffs->info + nInfoChars, line);
            nInfoChars += len;
        }
    };
    int assignedItems = sscanf(line, "%d %d %d %d %d", &minN, &maxN, &nTimes, &splineOrder, &splineSteps);
    if (assignedItems < 5)
    {
        fclose(f);
        return SHC_FILE_CONTENTS;
    }

	coeffs->minimumN = minN;
	coeffs->maximumN = maxN;
	coeffs->numberOfTimes = nTimes;
	coeffs->bSplineOrder = splineOrder;
	coeffs->bSplineSteps = splineSteps;

	nTerms = gsl_sf_legendre_array_n(maxN);
	coeffs->numberOfTerms = nTerms;

    nCoeffs = maxN * (maxN + 2) - (minN-1) * (minN - 1 +2);

	coeffs->polynomials = (double *)calloc(nTerms, sizeof(double));
	coeffs->derivatives = (double *)calloc(nTerms, sizeof(double));
	coeffs->aoverrpowers = (double *)calloc(maxN,  sizeof(double));
	coeffs->times = (double*)calloc(nTimes, sizeof(double));
	// Uses more memory than needed for coefficients, particularly in case of static field. 
    // Could be revised
	coeffs->gTimeSeries = (double*)calloc(nCoeffs * nTimes, sizeof(double));
	coeffs->hTimeSeries = (double*)calloc(nCoeffs * nTimes, sizeof(double));
	coeffs->gNow = (double*)calloc(nCoeffs, sizeof(double));
	coeffs->hNow = (double*)calloc(nCoeffs, sizeof(double));
	if (coeffs->polynomials == NULL || coeffs->aoverrpowers == NULL || coeffs->times == NULL || coeffs->gTimeSeries == NULL || coeffs->hTimeSeries == NULL || coeffs->gNow == NULL || coeffs->hNow == NULL)
	{
        status = SHC_MEMORY;
	}

	int il, im;
	for (int i = 0; i < nTimes; i++)
	{
		assignedItems = fscanf(f, "%lf", coeffs->times + i);
        if (assignedItems < 1)
        {
            fclose(f);
            return SHC_FILE_CONTENTS;
        }
	}
	size_t gRead = 0;
	size_t hRead = 0;
    size_t gCoeffs = 0;
    size_t hCoeffs = 0;
    int m = 0;
    size_t coeffCounter = 0;
    while (coeffCounter < nCoeffs)
    {
        assignedItems = fscanf(f, "%d %d", &il, &im);
        if (assignedItems < 2)
        {
            fclose(f);
            return SHC_FILE_CONTENTS;
        }

		if (im >= 0)
		{
			for (int i = 0; i < nTimes; i++)
			{
				assignedItems = fscanf(f, "%lf", coeffs->gTimeSeries + gRead);
				gRead++;
                if (assignedItems < 1)
                {
                    fclose(f);
                    return SHC_FILE_CONTENTS;
                }
			}
            gCoeffs++;
		}
		else
		{
			for (int i = 0; i < nTimes; i++)
			{
				assignedItems = fscanf(f, "%lf", coeffs->hTimeSeries + hRead);
				hRead++;
                if (assignedItems < 1)
                {
                    fclose(f);
                    return SHC_FILE_CONTENTS;
                }
			}
            hCoeffs++;
		}
        coeffCounter++;
	}

    if (coeffCounter != nCoeffs)
        return SHC_NUMBER_OF_COEFFICIENTS;

    coeffs->gCoeffs = gCoeffs;
    coeffs->hCoeffs = hCoeffs;

    fclose(f);

    coeffs->initialized = true;

	return SHC_OK;

}


void freeChaosCoefficients(ChaosCoefficients *coeffs)
{
    freeSHCCoefficients(&coeffs->core);
    freeSHCCoefficients(&coeffs->coreExtrapolation);
    freeSHCCoefficients(&coeffs->crust);

    return;
}

void freeSHCCoefficients(SHCCoefficients *coeffs)
{    
	if (coeffs->polynomials != NULL)
        free(coeffs->polynomials);
	if (coeffs->derivatives != NULL)
        free(coeffs->derivatives);
	if (coeffs->aoverrpowers != NULL)
        free(coeffs->aoverrpowers);
	if (coeffs->times != NULL)
        free(coeffs->times);
	if (coeffs->gTimeSeries != NULL)
        free(coeffs->gTimeSeries);
	if (coeffs->hTimeSeries != NULL)
        free(coeffs->hTimeSeries);
	if (coeffs->gNow != NULL)
        free(coeffs->gNow);
	if (coeffs->hNow != NULL)
        free(coeffs->hNow);

    return;
}

int interpolateSHCCoefficients(ChaosCoefficients *coeffs, double decimalYear)
{
    // Core coefficients at decimalYear (see decimalYearFromUnixTime), into core.gNow and core.hNow.
    // The core file samples each polynomial piece of the CHAOS B-spline at bSplineOrder times
    // (bSplineSteps apart), so Lagrange interpolation through the piece containing decimalYear
    // reproduces the spline. After the last core time, the extrapolation file (order 2, step 1)
    // gives the linear continuation. Before the first core time, the first snapshot is held.
    SHCCoefficients *c = &coeffs->core;
    double t = decimalYear;
    if (t > c->times[c->numberOfTimes - 1])
        c = &coeffs->coreExtrapolation;
    else if (t < c->times[0])
        t = c->times[0];

    int nTimes = c->numberOfTimes;
    int order = c->bSplineOrder;

    // First snapshot of the piece containing t; past the end, use the last piece
    int k = 0;
    while (k < nTimes - 1 && c->times[k + 1] <= t)
        k++;
    int first = (k / c->bSplineSteps) * c->bSplineSteps;
    if (first + order > nTimes)
        first = nTimes - order;

    double weights[SHC_MAX_SPLINE_ORDER];
    for (int p = 0; p < order; p++)
    {
        weights[p] = 1.0;
        for (int q = 0; q < order; q++)
        {
            if (q != p)
                weights[p] *= (t - c->times[first + q]) / (c->times[first + p] - c->times[first + q]);
        }
    }

    for (size_t i = 0; i < c->gCoeffs; i++)
    {
        coeffs->core.gNow[i] = 0.0;
        for (int p = 0; p < order; p++)
            coeffs->core.gNow[i] += weights[p] * c->gTimeSeries[i * nTimes + first + p];
    }
    for (size_t i = 0; i < c->hCoeffs; i++)
    {
        coeffs->core.hNow[i] = 0.0;
        for (int p = 0; p < order; p++)
            coeffs->core.hNow[i] += weights[p] * c->hTimeSeries[i * nTimes + first + p];
    }

    return SHC_OK;
}

// CHAOS decimal year, 2000 + MJD2000 / 365.25, the time scale of the SHC times and spline knots
double decimalYearFromUnixTime(double unixTime)
{
    return 2000.0 + (unixTime - 946684800.0) / 86400.0 / 365.25;
}
