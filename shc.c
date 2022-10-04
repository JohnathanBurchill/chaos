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

	coeffs->polynomials = (double *)calloc(nTerms, sizeof(double));
	coeffs->derivatives = (double *)calloc(nTerms, sizeof(double));
	coeffs->aoverrpowers = (double *)calloc(maxN,  sizeof(double));
	coeffs->times = (double*)calloc(nTimes, sizeof(double));
	// Uses more memory than needed for hnm. Could be revised
	coeffs->gTimeSeries = (double*)calloc(maxN * (maxN + 2) * nTimes, sizeof(double));
	coeffs->hTimeSeries = (double*)calloc(maxN * (maxN + 2) * nTimes, sizeof(double));
	coeffs->gNow = (double*)calloc(maxN * (maxN + 2), sizeof(double));
	coeffs->hNow = (double*)calloc(maxN * (maxN + 2), sizeof(double));
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
	for (int i = minN * (minN + 2); i <= maxN * (maxN+2); i++)
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
                if (assignedItems < 1)
                {
                    fclose(f);
                    return SHC_FILE_CONTENTS;
                }

				gRead++;
			}
		}
		else 
		{
			for (int i = 0; i < nTimes; i++)
			{
				assignedItems = fscanf(f, "%lf", coeffs->hTimeSeries + hRead);
                if (assignedItems < 1)
                {
                    fclose(f);
                    return SHC_FILE_CONTENTS;
                }
				hRead++;
			}
		}
	}

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

int interpolateSHCCoefficients(ChaosCoefficients *coeffs, int year, int month, int day)
{
	// Interpolate or extrapolate core model coefficients to the current day.
	// A resolution of 1 day is sufficient for CHAOS core model calculations for space physics

    double fractionalYear = 0.0;
	int status = 0;
    
    status = yearFraction(year, month, day, &fractionalYear);
	if (status != SHC_OK)
		return status;


	ssize_t coefficientTimeIndex = 0;
	ssize_t coefficientTimeIndexPlus1 = 0;
	double timeFraction = 0.0;
	double deltaTime = 0.0;

    int minN = 0;
    int maxN = 0;
    int nTimes = 0;
    double *times = NULL;
    double *gnmNow = coeffs->core.gNow;
    double *hnmNow = coeffs->core.hNow;
    double *gnm = NULL;
    double *hnm = NULL;
    // Linear extrapolation if time is greater than core field max time
    if (fractionalYear > coeffs->core.times[coeffs->core.numberOfTimes-1])
    {
        nTimes = coeffs->coreExtrapolation.numberOfTimes;
        coefficientTimeIndex = nTimes - 1;
        times = coeffs->coreExtrapolation.times;
        minN = coeffs->coreExtrapolation.minimumN;
        maxN = coeffs->coreExtrapolation.maximumN;
        gnm = coeffs->coreExtrapolation.gTimeSeries;
        hnm = coeffs->coreExtrapolation.hTimeSeries;

        // Only two times for core extrapolation coefficients
        if (nTimes != 2)
            return SHC_INTERPOLATION;

        // The first time is the same as the last time of the core coefficients
        deltaTime = times[1] - times[0];
        timeFraction = (fractionalYear - times[0]) / deltaTime;

        for (int i = (minN-1)*((minN-1) + 2); i < maxN * (maxN + 2); i++)
        {
            gnmNow[i] = gnm[i*nTimes + coefficientTimeIndex] + timeFraction * (gnm[i*nTimes + coefficientTimeIndexPlus1] - gnm[i*nTimes + coefficientTimeIndex]);
            // Calculates some irrelevant values at the end...
            hnmNow[i] = hnm[i*nTimes + coefficientTimeIndex] + timeFraction * (hnm[i*nTimes + coefficientTimeIndexPlus1] - hnm[i*nTimes + coefficientTimeIndex]);
        }

    }
    else
    {
        nTimes = coeffs->coreExtrapolation.numberOfTimes;
        coefficientTimeIndex = nTimes - 1;
        times = coeffs->coreExtrapolation.times;
        minN = coeffs->coreExtrapolation.minimumN;
        maxN = coeffs->coreExtrapolation.maximumN;
        gnm = coeffs->coreExtrapolation.gTimeSeries;
        hnm = coeffs->coreExtrapolation.hTimeSeries;

        // Linear interpolation
        while (coefficientTimeIndex > 0 && times[coefficientTimeIndex] > fractionalYear)
        {
            coefficientTimeIndex--;
        };

        // Linear interpolation is sufficient given the number of model times is 5x the knot points.
        if ((coefficientTimeIndex == (nTimes - 1)) || (times[coefficientTimeIndex] > fractionalYear))
        {
            coefficientTimeIndexPlus1 = coefficientTimeIndex;
            deltaTime = 1.0; // arbitrary, as we divide this into 0.0
        }
        else
        {
            coefficientTimeIndexPlus1 = coefficientTimeIndexPlus1 + 1;
            deltaTime = times[coefficientTimeIndexPlus1] - times[coefficientTimeIndex];
        }
        timeFraction = (fractionalYear - times[coefficientTimeIndex]) / deltaTime;
        for (int i = (minN-1)*((minN-1) + 2); i < maxN * (maxN + 2); i++)
        {
            gnmNow[i] = gnm[i*nTimes + coefficientTimeIndex] + timeFraction * (gnm[i*nTimes + coefficientTimeIndexPlus1] - gnm[i*nTimes + coefficientTimeIndex]);
            // Calculates some irrelevant values at the end...
            hnmNow[i] = hnm[i*nTimes + coefficientTimeIndex] + timeFraction * (hnm[i*nTimes + coefficientTimeIndexPlus1] - hnm[i*nTimes + coefficientTimeIndex]);
            // printf("%lf %lf\n", gnmNow[i], hnmNow[i]);
        }


    }

	// Crustal field is static, so copy g and h into gNow and h Now
	// simply loops over all indices, but those for n < 21 are irrelevant.
	// If there is more than 1 time for the crustal field, abort
    nTimes = coeffs->crust.numberOfTimes;
	if (nTimes != 1)
	{
		fprintf(stderr, "Expected 1 crustal field time, got %d times\n", nTimes);
		return SHC_INTERPOLATION;
	}
	for (int i = 0; i < coeffs->crust.numberOfTerms; i++)
	{
		coeffs->crust.gNow[i] = coeffs->crust.gTimeSeries[i];
		coeffs->crust.hNow[i] = coeffs->crust.hTimeSeries[i];
	}

    return SHC_OK;
}

// Calculates day of year: 1 January is day 1.
int yearFraction(long year, long month, long day, double* fractionalYear)
{
    time_t date;
    struct tm dateStruct;
    dateStruct.tm_year = year - 1900;
    dateStruct.tm_mon = month - 1;
    dateStruct.tm_mday = day;
    dateStruct.tm_hour = 0;
    dateStruct.tm_min = 0;
    dateStruct.tm_sec = 0;
    dateStruct.tm_yday = 0;
    date = timegm(&dateStruct);
    struct tm *dateStructUpdated = gmtime(&date);
    if (dateStructUpdated == NULL)
        return SHC_FRACTIONAL_YEAR;
	// How does this work for leap years?
    *fractionalYear = (double) year + (double)(dateStructUpdated->tm_yday + 1)/365.25;
    return SHC_OK;
}
