/*

    CHAOS: chaos.c

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

#include "cdf_utils.h"
#include "chaos_settings.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <signal.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_math.h>

#define NMAGVARS 5

#define MAG_INPUT_DATASET "LR_1B"
#define INTERPOLATION_SKIP 4 // 4 s
//#define MAG_INPUT_DATASET "HR_1B"
//#define INTERPOLATION_SKIP 200 // 4 s

// declare functions

// Handle Ctrl-C
// From https://stackoverflow.com/questions/4217037/catch-ctrl-c-in-c
volatile sig_atomic_t keep_running = 1;
static void sig_handler(int ignored)
{
    (void)ignored; // Avoid warning
    keep_running = 0;
}

void usage(const char *name);
int yearFraction(long year, long month, long day, double* fractionalYear);
int calculateField(double r, double theta, double phi, int minN, int maxN, double *gnm, double *hnm, double *polynomials, double *derivatives, double *aoverrpowers, double *bn, double *be, double *bc);

int loadModel(const char *filename, int *minimumN, int *maximumN, int *numberOfTimes, size_t *numberOfTerms, int *bSplineOrder, int *bSplineSteps, double **times, double **gTimeSeries, double **hTimeSeries, double **gNow, double **hNow, double **polynomials, double **derivatives, double **aoverrpowers);


char infoHeader[50];

int main (int argc, char **argv)
{

	// Handle Ctrl-C
	signal(SIGINT, sig_handler);

	time_t processingStartTime = time(NULL);

	int status = 0;

	printf("Hi\n");

	double degrees = M_PI / 180.0;
	double a = 6371.2;
	double inputTime = 0.;
	double r = 0.;
	double theta = 0.0 * degrees;
	double phi = 0.0 * degrees;

	int nTimes = 0;
	int maxN = 0;
	int minN = 0;
	int bSplineOrder = 0;
	int bSplineSteps = 0;
	size_t nTerms = 0;
	double *polynomials = NULL;
	double *derivatives = NULL;
	double *aoverrpowers = NULL;
	double *times = NULL;
	double *gnm = NULL;
	double *hnm = NULL;
	double *gnmNow = NULL;
	double *hnmNow = NULL;

	int nTimesCrustal = 0;
	int maxNCrustal = 0;
	int minNCrustal = 0;
	int bSplineOrderCrustal = 0;
	int bSplineStepsCrustal = 0;
	size_t nTermsCrustal = 0;
	double *polynomialsCrustal = NULL;
	double *derivativesCrustal = NULL;
	double *aoverrpowersCrustal = NULL;
	double *timesCrustal = NULL;
	double *gnmCrustal = NULL;
	double *hnmCrustal = NULL;
	double *gnmNowCrustal = NULL;
	double *hnmNowCrustal = NULL;


	double *bCore = NULL;
	double *bCrust = NULL;
	double *dbMeas = NULL;

	size_t nInputs = 0;
	uint8_t *magVariables[NMAGVARS];

	// Timestamp, Latitude, Longitude, Radius, B_NEC
	for (int i = 0; i < NMAGVARS; i++)
	{
		magVariables[i] = NULL;
	}

	char coreFile[FILENAME_MAX];
	char crustalFile[FILENAME_MAX];
	char magFilename[FILENAME_MAX];

	char outputFilename[FILENAME_MAX];

	if (argc != 5)
	{
		usage(argv[0]);
		goto cleanup;
	}
	char *satDate = argv[1];
	char *coeffDir = argv[2];
	char *magDir = argv[3];
	char *outputDir = argv[4];

	char satellite = satDate[0];
	long year, month, day;
    int valuesRead = sscanf(satDate+1, "%4ld%2ld%2ld", &year, &month, &day);
    if (valuesRead != 3)
    {
        fprintf(stdout, "CHAOS called as:\n \"%s %s %s %s %s\"\n Unable to parse date \"\". Exiting.\n", argv[0], argv[1], argv[2], argv[3], argv[4]);
		goto cleanup;
    }
	sprintf(infoHeader, "CHAOS %s: ", satDate);

	double beginTime;
	double endTime;
	status = getOutputFilename(satellite, year, month, day, outputDir, &beginTime, &endTime, outputFilename);
	if (status != 0)
	{
		printf("Could not get output filename.\n");
		goto cleanup;
	}
	char fullOutputFilename[FILENAME_MAX];
	sprintf(fullOutputFilename, "%s", outputFilename);
	sprintf(fullOutputFilename+strlen(outputFilename), "%s", ".cdf");
	if (access(fullOutputFilename, F_OK) == 0)
	{
		printf("%sOutput CDF file exists. Exiting.\n", infoHeader);
		goto cleanup;
	}

	// Core model coefficients
	sprintf(coreFile, "%s/%s", coeffDir, "CHAOS-7.9_core.shc");
	printf("Core model file: %s\n", coreFile);
	status = loadModel(coreFile, &minN, &maxN, &nTimes, &nTerms, &bSplineOrder, &bSplineSteps, &times, &gnm, &hnm, &gnmNow, &hnmNow, &polynomials, &derivatives, &aoverrpowers);
	if (status != 0)
	{
		printf("Could not load core field model coefficients.\n");
		goto cleanup;
	}

	// Interpolate core model coefficients to the current day.
	// A resolution of 1 day will be sufficient for CHAOS core model calculations
	// Given that inputs are from daily MAG_HR CDF files, we only need
	// to interpolate once.
	double fractionalYear = 0.0;
	status = yearFraction(year, month, day, &fractionalYear);
	if (status != 0)
	{
		printf("Could not calculate fractional year.\n");
		goto cleanup;
	}
	printf("Factional year: %.4lf\n", fractionalYear);
	size_t coefficientTimeIndex = nTimes - 1;
	size_t coefficientTimeIndexPlus1 = 0;
	double timeFraction = 0.0;
	double deltaTime = 0.0;
	while (times[coefficientTimeIndex] > fractionalYear && coefficientTimeIndex > 0)
	{
		coefficientTimeIndex--;
	};
	// Linear interpolation is sufficient given the number of model times is 5x the knot points.
	// Constant extrapolation, should be revised to linearly extrapolate using dedicated model coefficients.
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
	printf ("maxN: %d\n", maxN);
	for (int i = (minN-1)*((minN-1) + 2); i < maxN * (maxN + 2); i++)
	{
		gnmNow[i] = gnm[i*nTimes + coefficientTimeIndex] + timeFraction * (gnm[i*nTimes + coefficientTimeIndexPlus1] - gnm[i*nTimes + coefficientTimeIndex]);
		// Calculates some irrelevant values at the end...
		hnmNow[i] = hnm[i*nTimes + coefficientTimeIndex] + timeFraction * (hnm[i*nTimes + coefficientTimeIndexPlus1] - hnm[i*nTimes + coefficientTimeIndex]);
		// printf("%lf %lf\n", gnmNow[i], hnmNow[i]);
	}

	// Crustal model coefficients
	sprintf(crustalFile, "%s/%s", coeffDir, "CHAOS-7.9_static.shc");
	printf("Crustal model file: %s\n", crustalFile);
	status = loadModel(crustalFile, &minNCrustal, &maxNCrustal, &nTimesCrustal, &nTermsCrustal, &bSplineOrderCrustal, &bSplineStepsCrustal, &timesCrustal, &gnmCrustal, &hnmCrustal, &gnmNowCrustal, &hnmNowCrustal, &polynomialsCrustal, &derivativesCrustal, &aoverrpowersCrustal);
	if (status != 0)
	{
		printf("Could not load crustal field model coefficients.\n");
		goto cleanup;
	}
	// Crustal field is static, so copy g and h into gNow and h Now
	// simply loops over all indices, but those for n < 21 are irrelevant.
	// If there is more than 1 time for the crustal field, abort
	if (nTimesCrustal != 1)
	{
		printf("Expected 1 crustal field time, got %d\n", nTimesCrustal);
		goto cleanup;
	}
	for (int i = 0; i < nTermsCrustal; i++)
	{
		gnmNowCrustal[i] = gnmCrustal[i];
		hnmNowCrustal[i] = hnmCrustal[i];
	}

	char *magVariableNames[NMAGVARS] = {
		"Timestamp",
		"Latitude",
		"Longitude",
		"Radius",
		"B_NEC"
	};
	// Magnetic field input data
	// LR_1B product for development, much faster load time than HR_1B
	if (getInputFilename(satellite, year, month, day, magDir, MAG_INPUT_DATASET, magFilename))
	// if (getInputFilename(satellite, year, month, day, magDir, "HR_1B", magFilename))
    {
        fprintf(stdout, "%sMAG HR_1B input file is not available. Exiting.\n", infoHeader);
        exit(1);
    }
	printf("%sReading inputs from %s\n", infoHeader, magFilename);
	loadCdf(magFilename, magVariableNames, NMAGVARS, magVariables, &nInputs);
	if (nInputs == 0)
	{
		printf("Found no inputs from MAG_HR file.\n");
		goto cleanup;
	}

	// Measured fields
	dbMeas = (double*)malloc(nInputs * 3 * sizeof(double));

	// Model fields
	bCore = (double*)malloc(nInputs * 3 * sizeof(double));
	bCrust = (double*)malloc(nInputs * 3 * sizeof(double));
	if (dbMeas == NULL || bCore == NULL || bCrust == NULL)
	{
		printf("Memory issue.\n");
		goto cleanup;
	}

	double deltaT = 0.0;
	double interpolationFraction = 1.0;
	size_t lastIndex = 0;

	if (keep_running == 1)
	{
		printf("Calculating fields...\n");


		// Get fields for first time
		inputTime = ((double*)magVariables[0])[0];
		theta = (90.0 - ((double*)magVariables[1])[0]) * degrees;
		phi = ((double*)magVariables[2])[0] * degrees;
		r = ((double*)magVariables[3])[0]/1000.;

		status = calculateField(r, theta, phi, minN, maxN, gnmNow, hnmNow, polynomials, derivatives, aoverrpowers, bCore, bCore+1, bCore+2);
		if (status != 0)
		{
			printf("Could not calculate core field\n");
			goto cleanup;
		}
		status = calculateField(r, theta, phi, minNCrustal, maxNCrustal, gnmNowCrustal, hnmNowCrustal, polynomialsCrustal, derivativesCrustal, aoverrpowersCrustal, bCrust, bCrust+1, bCrust+2);
		if (status != 0)
		{
			printf("Could not calculate crustal field\n");
			goto cleanup;
		}
		dbMeas[0] = ((double*)magVariables[4])[0] - bCore[0] - bCrust[0];
		dbMeas[1] = ((double*)magVariables[4])[1] - bCore[1] - bCrust[1];
		dbMeas[2] = ((double*)magVariables[4])[2] - bCore[2] - bCrust[2];

		lastIndex = 0;

		for (size_t t = INTERPOLATION_SKIP; t < nInputs && keep_running == 1; t += INTERPOLATION_SKIP)
		// for (size_t t = 0; t < nInputs && keep_running == 1; t+=1*60*5)
		// for (size_t t = 0; t < nInputs && keep_running == 1; t+=50*60*5)
		{
			inputTime = ((double*)magVariables[0])[t];
			theta = (90.0 - ((double*)magVariables[1])[t]) * degrees;
			phi = ((double*)magVariables[2])[t] * degrees;
			r = ((double*)magVariables[3])[t]/1000.;

			status = calculateField(r, theta, phi, minN, maxN, gnmNow, hnmNow, polynomials, derivatives, aoverrpowers, bCore+t*3, bCore+t*3+1, bCore+t*3+2);
			if (status != 0)
			{
				printf("Could not calculate core field\n");
				goto cleanup;
			}

			status = calculateField(r, theta, phi, minNCrustal, maxNCrustal, gnmNowCrustal, hnmNowCrustal, polynomialsCrustal, derivativesCrustal, aoverrpowersCrustal, bCrust+t*3, bCrust+t*3+1, bCrust+t*3+2);
			if (status != 0)
			{
				printf("Could not calculate crustal field\n");
				goto cleanup;
			}
			dbMeas[t*3]   = ((double*)magVariables[4])[t*3 + 0] - bCore[t*3 + 0] - bCrust[t*3 + 0];
			dbMeas[t*3+1] = ((double*)magVariables[4])[t*3 + 1] - bCore[t*3 + 1] - bCrust[t*3 + 1];
			dbMeas[t*3+2] = ((double*)magVariables[4])[t*3 + 2] - bCore[t*3 + 2] - bCrust[t*3 + 2];

			// Linearly interpolate at skipped epochs
			deltaT = ((double*)magVariables[0])[t] - ((double*)magVariables[0])[t - INTERPOLATION_SKIP];
			if (deltaT <= 0.)
				deltaT = 1.0; // arbitrary
			for (size_t i = t-INTERPOLATION_SKIP + 1; i < t; i++)
			{
				interpolationFraction = (((double*)magVariables[0])[i] - ((double*)magVariables[0])[t-INTERPOLATION_SKIP]) / deltaT;
				bCore[i*3] = bCore[(t-INTERPOLATION_SKIP)*3] + interpolationFraction * (bCore[t*3] - bCore[(t-INTERPOLATION_SKIP)*3]);
				bCore[i*3+1] = bCore[(t-INTERPOLATION_SKIP)*3+1] + interpolationFraction * (bCore[t*3+1] - bCore[(t-INTERPOLATION_SKIP)*3+1]);
				bCore[i*3+2] = bCore[(t-INTERPOLATION_SKIP)*3+2] + interpolationFraction * (bCore[t*3+2] - bCore[(t-INTERPOLATION_SKIP)*3+2]);

				bCrust[i*3] = bCrust[(t-INTERPOLATION_SKIP)*3] + interpolationFraction * (bCrust[t*3] - bCrust[(t-INTERPOLATION_SKIP)*3]);
				bCrust[i*3+1] = bCrust[(t-INTERPOLATION_SKIP)*3+1] + interpolationFraction * (bCrust[t*3+1] - bCrust[(t-INTERPOLATION_SKIP)*3+1]);
				bCrust[i*3+2] = bCrust[(t-INTERPOLATION_SKIP)*3+2] + interpolationFraction * (bCrust[t*3+2] - bCrust[(t-INTERPOLATION_SKIP)*3+2]);

				// Copy skipped measured fields
				dbMeas[i*3]   = ((double*)magVariables[4])[i*3 + 0] - bCore[i*3 + 0] - bCrust[i*3 + 0];
				dbMeas[i*3+1] = ((double*)magVariables[4])[i*3 + 1] - bCore[i*3 + 1] - bCrust[i*3 + 1];
				dbMeas[i*3+2] = ((double*)magVariables[4])[i*3 + 2] - bCore[i*3 + 2] - bCrust[i*3 + 2];

			}
			lastIndex = t;
			// printf("t=%ld\n", t);
			// With magnetic potential
			// printf("magneticPotential(time=%.1lf, r=%6.1lf, colatitude=%5.1lf, longitude=%5.1lf) = %.1lf, br = %.1lf, btheta = %.1lf nT\n", time, r, theta/degrees, phi/degrees, magneticPotential, br, btheta);
			// Field vectors only
	//		printf("time=%.1lf, r=%6.1lf, latitude=%5.1lf, longitude=%6.1lf: (%8.1lf, %8.1lf, %8.1lf) nT (NEC)\n", inputTime, r, 90.0 - theta/degrees, phi/degrees, -btheta, bphi, -br);

			// printf("time=%.1lf: model/measured=(%8.1lf/%8.1lf, %8.1lf/%8.1lf, %8.1lf/%8.1lf) nT (NEC)\n", inputTime, -btheta, bn, bphi, be, -br, bc);
			// Delta-B (core removed)
			// printf("time=%.1lf: DeltaB = (%8.1lf, %8.1lf, %8.1lf) nT (NEC)\n", inputTime, bn - bcn, be - bce, bc - bcc);
			// Delta-B (core and crust removed)
			// printf("time=%.1lf: DeltaB = (%8.1lf, %8.1lf, %8.1lf) nT (NEC)\n", inputTime, bn - bcn - bcrn, be - bce - bcre, bc - bcc - bcrc);
			// B crustal
			// printf("time=%.1lf: Crustal B = (%8.1lf, %8.1lf, %8.1lf) nT (NEC)\n", inputTime, bcrn, bcre, bcrc);
			// Max crustal magnitudes

		}

		// Calculate actual field for remaining inputs (up to INTEGRATION_SKIP-1 of them)
		for (size_t t = lastIndex+1; t < nInputs && keep_running == 1; t ++)
		{
			inputTime = ((double*)magVariables[0])[t];
			theta = (90.0 - ((double*)magVariables[1])[t]) * degrees;
			phi = ((double*)magVariables[2])[t] * degrees;
			r = ((double*)magVariables[3])[t]/1000.;

			status = calculateField(r, theta, phi, minN, maxN, gnmNow, hnmNow, polynomials, derivatives, aoverrpowers, bCore+t*3, bCore+t*3+1, bCore+t*3+2);
			if (status != 0)
			{
				printf("Could not calculate core field\n");
				goto cleanup;
			}

			status = calculateField(r, theta, phi, minNCrustal, maxNCrustal, gnmNowCrustal, hnmNowCrustal, polynomialsCrustal, derivativesCrustal, aoverrpowersCrustal, bCrust+t*3, bCrust+t*3+1, bCrust+t*3+2);
			if (status != 0)
			{
				printf("Could not calculate crustal field\n");
				goto cleanup;
			}

			dbMeas[t*3]   = ((double*)magVariables[4])[t*3 + 0] - bCore[t*3 + 0] - bCrust[t*3 + 0];
			dbMeas[t*3+1] = ((double*)magVariables[4])[t*3 + 1] - bCore[t*3 + 1] - bCrust[t*3 + 1];
			dbMeas[t*3+2] = ((double*)magVariables[4])[t*3 + 2] - bCore[t*3 + 2] - bCrust[t*3 + 2];

		}

	}
	if (keep_running == 0)
	{
		printf("%sInterrupted (SIGINT).\n", infoHeader);
		goto cleanup;
	}

	// printf("Max Crustal B magnitudes = (%8.1lf, %8.1lf, %8.1lf) nT (NEC)\n", maxBcrn, maxBcre, maxBcrc);

	// Export CDF of times, measured B, core B, and crustal B
	status = exportCdf(outputFilename, satellite, EXPORT_VERSION_STRING, (double*)magVariables[0], bCore, bCrust, dbMeas, nInputs);
	if (status != 0)
	{
		printf("Could not export fields.\n");
	}

	time_t processingStopTime = time(NULL);
	exportMetaInfo(outputFilename, magFilename, coreFile, crustalFile, nInputs, processingStartTime, processingStopTime);

cleanup:
	if (polynomials != NULL) free(polynomials);
	if (derivatives != NULL) free(derivatives);
	if (aoverrpowers != NULL) free(aoverrpowers);
	if (times != NULL) free(times);
	if (gnm != NULL) free(gnm);
	if (hnm != NULL) free(hnm);
	if (gnmNow != NULL) free(gnmNow);
	if (hnmNow != NULL) free(hnmNow);
	if (polynomialsCrustal != NULL) free(polynomialsCrustal);
	if (derivativesCrustal != NULL) free(derivativesCrustal);
	if (aoverrpowersCrustal != NULL) free(aoverrpowersCrustal);
	if (timesCrustal != NULL) free(timesCrustal);
	if (gnmCrustal != NULL) free(gnmCrustal);
	if (hnmCrustal != NULL) free(hnmCrustal);
	if (gnmNowCrustal != NULL) free(gnmNowCrustal);
	if (hnmNowCrustal != NULL) free(hnmNowCrustal);
	for (int i = 0; i < NMAGVARS; i++)
	{
		if (magVariables[i] != NULL) free(magVariables[i]);
	}
	if (dbMeas != NULL) free(dbMeas);
	if (bCore != NULL) free(bCore);
	if (bCrust != NULL) free(bCrust);

	return 0;
}

void usage(const char* name)
{
	printf("Usage: %s XYYYYMMDD chaosModelCoefficientsDir magCdfDir outputDir\n", name);
	printf(" X: satellite letter A, B, or C\n");
	printf(" YYYYMMDD: year, month, day\n");
	printf(" chaosModelCoefficientsDir: directory containing SHC files\n");
	printf(" magCdfDir: directory containing MAG_HR CDFs\n");
	printf(" outputDir: directory to store magnetic field vectors\n");
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
    {
        fprintf(stdout, "%sUnable to get day of year from specified date.\n", infoHeader);
        *fractionalYear = 0;
        return -1;
    }
	// How does this work for leap years?
    *fractionalYear = (double) year + (double)(dateStructUpdated->tm_yday + 1)/365.25;
    return 0;
}

int calculateField(double r, double theta, double phi, int minN, int maxN, double *gnm, double *hnm, double *polynomials, double *derivatives, double *aoverrpowers, double *bn, double *be, double *bc)
{
	double a = 6371.2;
	double aoverr = a/r;
	double magneticPotential = 0.0;
	double magneticPotentialN = 0.0;
	double magneticPotentialForThetaN = 0.0;
	double magneticPotentialForPhiN = 0.0;
	double br = 0.0;
	double btheta = 0.0;
	double bphi = 0.0;
	size_t gRead = 0;
	size_t hRead = 0;

	int status = 0;

	aoverrpowers[0] = aoverr * aoverr; // (a/r)^n+1, n starting at 1
	for (int i = 2; i <= maxN; i++)
	{
		aoverrpowers[i-1] = aoverrpowers[i-2] * aoverr;
	}

	status = gsl_sf_legendre_deriv_array(GSL_SF_LEGENDRE_SCHMIDT, maxN, cos(theta), polynomials, derivatives);
	if (status)
	{
		printf("status: %s\n", gsl_strerror(status));
		return -1;
	}

	for (int n = minN; n <= maxN; n++)
	{
		magneticPotentialN = gnm[gRead] * polynomials[gsl_sf_legendre_array_index(n, 0)];
		magneticPotentialForThetaN = gnm[gRead] * derivatives[gsl_sf_legendre_array_index(n, 0)] * (-sin(theta));
		magneticPotentialForPhiN = 0.0;
		gRead++;
		for (int m = 1; m <= n; m++)
		{
			magneticPotentialN += (gnm[gRead] * cos((double)m*phi) + hnm[hRead] * sin((double)m*phi) ) * polynomials[gsl_sf_legendre_array_index(n, m)];
			magneticPotentialForThetaN += (gnm[gRead] * cos((double)m*phi) + hnm[hRead] * sin((double)m*phi) ) * derivatives[gsl_sf_legendre_array_index(n, m)] * (-sin(theta));
			magneticPotentialForPhiN += (gnm[gRead] * (-m*sin((double)m*phi)) + hnm[hRead] * m*cos((double)m*phi) ) * polynomials[gsl_sf_legendre_array_index(n, m)];
			gRead++;
			hRead++;
		}
		magneticPotentialN *= aoverrpowers[n-1] * a;
		magneticPotentialForThetaN *= aoverrpowers[n-1] * a;
		magneticPotentialForPhiN *= aoverrpowers[n-1] * a;

		magneticPotential += magneticPotentialN;
		// Br = -di potential by di r
		br += -magneticPotentialN * (-(n+1.) / r);
		// Btheta = -di potential by di theta / r
		btheta += -magneticPotentialForThetaN / r;
		// Bphi = -di potential by di phi / (r*sin(theta))
		bphi += -magneticPotentialForPhiN / (r * sin(theta));
	}

	*bn = -btheta;
	*be = bphi;
	*bc = -br;

	return 0;

}

int loadModel(const char *filename, int *minimumN, int *maximumN, int *numberOfTimes, size_t *numberOfTerms, int *bSplineOrder, int *bSplineSteps, double **times, double **gTimeSeries, double **hTimeSeries, double **gNow, double **hNow, double **polynomials, double **derivatives, double **aoverrpowers)
{
	int status = 0;
	char line[255];

	FILE *f = fopen(filename, "r");	
	if (f == NULL)
	{
		printf("Could not open core SHC file.\n");
		status = -1;
		goto cleanup;
	}

	int minN = 0, maxN = 0, nTimes = 0, splineOrder = 0, splineSteps = 0; 
	size_t nTerms = 0;
	while (fgets(line, 255, f) != NULL && line[0] == '#');
	sscanf(line, "%d %d %d %d %d", &minN, &maxN, &nTimes, &splineOrder, &splineSteps);

	*minimumN = minN;
	*maximumN = maxN;
	*numberOfTimes = nTimes;
	*bSplineOrder = splineOrder;
	*bSplineSteps = splineSteps;

	nTerms = gsl_sf_legendre_array_n(maxN);
	*numberOfTerms = nTerms;

	*polynomials = (double *)malloc(nTerms * sizeof(double));
	*derivatives = (double *)malloc(nTerms * sizeof(double));
	*aoverrpowers = (double *)malloc(maxN * sizeof(double));
	*times = (double*)malloc(nTimes * sizeof(double));
	// Uses more memory than needed for hnm. Could be revised
	*gTimeSeries = (double*)malloc(maxN * (maxN + 2) * nTimes * sizeof(double));
	*hTimeSeries = (double*)malloc(maxN * (maxN + 2) * nTimes * sizeof(double));
	*gNow = (double*)malloc(maxN * (maxN + 2) * sizeof(double));
	*hNow = (double*)malloc(maxN * (maxN + 2) * sizeof(double));
	if (*polynomials == NULL || *aoverrpowers == NULL || *times == NULL || *gTimeSeries == NULL || *hTimeSeries == NULL || *gNow == NULL || *hNow == NULL)
	{
		printf("Cannot remember anything. Check my memory.\n");
		goto cleanup;
	}

	int il, im;
	int assignedItems = 0;
	for (int i = 0; i < nTimes; i++)
	{
		assignedItems = fscanf(f, "%lf", *times + i);
	}
	size_t gRead = 0;
	size_t hRead = 0;
	for (int i = minN * (minN + 2); i <= maxN * (maxN+2); i++)
	{
		assignedItems = fscanf(f, "%d %d", &il, &im);
		if (im >= 0)
		{
			for (int i = 0; i < nTimes; i++)
			{
				assignedItems = fscanf(f, "%lf", *gTimeSeries + gRead);
				gRead++;
			}
		}
		else 
		{
			for (int i = 0; i < nTimes; i++)
			{
				assignedItems = fscanf(f, "%lf", *hTimeSeries + hRead);
				hRead++;
			}
		}
	}
	// for (int i = 0; i < nTimes; i++)
	// {
	// 	printf("t = %lf\n", times[i]);
	// }
	// gRead = 0;
	// hRead = 0;
	// for (int l = minN; l <= maxN; l++)
	// {
	// 	for (int m = 0; m <= l; m++)
	// 	{
	// 		printf("g_%d_%d = %lf\n", l, m, gnm[gRead*nTimes + 0]);
	// 		gRead++;
	// 		if (m > 0)
	// 		{
	// 			printf("h_%d_%d = %lf\n", l, m, hnm[hRead*nTimes + 0]);
	// 			hRead++;
	// 		}
	// 	}
	// }

cleanup:	
	if (f != NULL) fclose(f);

	return 0;

}