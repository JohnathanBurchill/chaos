#include "cdf_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <signal.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_math.h>

#define NMAGVARS 5

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
int calculateField(double r, double theta, double phi, int maxN, double *gnm, double *hnm, double *polynomials, double *derivatives, double *aoverrpowers, double *bn, double *be, double *bc);


char infoHeader[50];

int main (int argc, char **argv)
{

	// Handle Ctrl-C
	signal(SIGINT, sig_handler);

	int status = 0;

	printf("Hi\n");

	double degrees = M_PI / 180.0;
	double a = 6371.2;
	double time = 0.;
	double r = 0.;
	double theta = 0.0 * degrees;
	double phi = 0.0 * degrees;

	double *polynomials = NULL;
	double *derivatives = NULL;
	double *aoverrpowers = NULL;
	double *times = NULL;
	double *gnm = NULL;
	double *hnm = NULL;
	double *gnmNow = NULL;
	double *hnmNow = NULL;

	size_t nInputs = 0;
	uint8_t *magVariables[NMAGVARS];

	// Timestamp, Latitude, Longitude, Radius, B_NEC
	for (int i = 0; i < NMAGVARS; i++)
	{
		magVariables[i] = NULL;
	}

	FILE *f = NULL;

	char coreFile[FILENAME_MAX];
	char magFilename[FILENAME_MAX];

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

	sprintf(coreFile, "%s/%s", coeffDir, "CHAOS-7.9_core.shc");
	printf("Core model file: %s\n", coreFile);
	char line[255];

	f = fopen(coreFile, "r");	
	if (f == NULL)
	{
		printf("Could not open core SHC file.\n");
		goto cleanup;
	}

	int minN, maxN, nTimes, bSplineOrder, bSplineSteps; 
	while (fgets(line, 255, f) != NULL && line[0] == '#');
	sscanf(line, "%d %d %d %d %d", &minN, &maxN, &nTimes, &bSplineOrder, &bSplineSteps);

	printf("MaxN: %d, nTimes: %d\n", maxN, nTimes);

	size_t nTerms = gsl_sf_legendre_array_n(maxN);
	polynomials = malloc(nTerms * sizeof(double));
	derivatives = malloc(nTerms * sizeof(double));
	aoverrpowers = malloc(maxN * sizeof(double));
	times = (double*)malloc(nTimes * sizeof(double));
	// Uses more memory than needed for hnm. Could be revised
	gnm = (double*)malloc(maxN * (maxN + 2) * nTimes * sizeof(double));
	hnm = (double*)malloc(maxN * (maxN + 2) * nTimes * sizeof(double));
	gnmNow = (double*)malloc(maxN * (maxN + 2) * sizeof(double));
	hnmNow = (double*)malloc(maxN * (maxN + 2) * sizeof(double));
	if (polynomials == NULL || aoverrpowers == NULL || times == NULL || gnm == NULL || hnm == NULL || gnmNow == NULL || hnmNow == NULL)
	{
		printf("Cannot remember anything. Check my memory.\n");
		goto cleanup;
	}

	int il, im;
	for (int i = 0; i < nTimes; i++)
	{
		fscanf(f, "%lf", times + i);
	}
	size_t gRead = 0;
	size_t hRead = 0;
	for (int i = 0; i < maxN * (maxN+2); i++)
	{
		fscanf(f, "%d %d", &il, &im);
		if (im >= 0)
		{
			for (int i = 0; i < nTimes; i++)
			{
				fscanf(f, "%lf", gnm + gRead);
				gRead++;
			}
		}
		else 
		{
			for (int i = 0; i < nTimes; i++)
			{
				fscanf(f, "%lf", hnm + hRead);
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
	// for (int l = 1; l <= maxN; l++)
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
	fclose(f);
	f = NULL;

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
	for (int i = 0; i < maxN * (maxN + 2); i++)
	{
		gnmNow[i] = gnm[i*nTimes + coefficientTimeIndex] + timeFraction * (gnm[i*nTimes + coefficientTimeIndexPlus1] - gnm[i*nTimes + coefficientTimeIndex]);
		// Calculates some irrelevant values at the end...
		hnmNow[i] = hnm[i*nTimes + coefficientTimeIndex] + timeFraction * (hnm[i*nTimes + coefficientTimeIndexPlus1] - hnm[i*nTimes + coefficientTimeIndex]);
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
	if (getInputFilename(satellite, year, month, day, magDir, "LR_1B", magFilename))
	// if (getInputFilename(satellite, year, month, day, magDir, "HR_1B", magFilename))
    {
        fprintf(stdout, "%sMAG HR_1B input file is not available. Exiting.\n", infoHeader);
        exit(1);
    }

	loadCdf(magFilename, magVariableNames, NMAGVARS, magVariables, &nInputs);
	if (nInputs == 0)
	{
		printf("Found no inputs from MAG_HR file.\n");
		goto cleanup;
	}

	// Measured fields
	double bn = 0.0;
	double be = 0.0;
	double bc = 0.0;

	// Core field
	double bcn = 0.0;
	double bce = 0.0;
	double bcc = 0.0;

	double aoverr = 1.0;

	if (keep_running == 1)
		printf("Calculating magnetic potential...\n");
	else
		printf("Calculation interrupted.\n");
	
	for (size_t t = 0; t < nInputs && keep_running == 1; t++)
	// for (size_t t = 0; t < nInputs && keep_running == 1; t+=1*60*5)
	// for (size_t t = 0; t < nInputs && keep_running == 1; t+=50*60*5)
	{
		time = ((double*)magVariables[0])[t];
	 	theta = (90.0 - ((double*)magVariables[1])[t]) * degrees;
		phi = ((double*)magVariables[2])[t] * degrees;
		r = ((double*)magVariables[3])[t]/1000.;

		status = calculateField(r, theta, phi, maxN, gnmNow, hnmNow, polynomials, derivatives, aoverrpowers, &bcn, &bce, &bcc);
		if (status != 0)
		{
			printf("Could not calculate core field\n");
			goto cleanup;
		}

		// printf("t=%ld\n", t);
		// With magnetic potential
		// printf("magneticPotential(time=%.1lf, r=%6.1lf, colatitude=%5.1lf, longitude=%5.1lf) = %.1lf, br = %.1lf, btheta = %.1lf nT\n", time, r, theta/degrees, phi/degrees, magneticPotential, br, btheta);
		// Field vectors only
//		printf("time=%.1lf, r=%6.1lf, latitude=%5.1lf, longitude=%6.1lf: (%8.1lf, %8.1lf, %8.1lf) nT (NEC)\n", time, r, 90.0 - theta/degrees, phi/degrees, -btheta, bphi, -br);
		// Compare model with measured fields
		bn = ((double*)magVariables[4])[t*3 + 0];
		be = ((double*)magVariables[4])[t*3 + 1];
		bc = ((double*)magVariables[4])[t*3 + 2];
		// printf("time=%.1lf: model/measured=(%8.1lf/%8.1lf, %8.1lf/%8.1lf, %8.1lf/%8.1lf) nT (NEC)\n", time, -btheta, bn, bphi, be, -br, bc);
		// Delta-B (core removed)
		printf("time=%.1lf: DeltaB = (%8.1lf, %8.1lf, %8.1lf) nT (NEC)\n", time, bn - bcn, be - bce, bc - bcc);

	}

cleanup:
	if (polynomials != NULL) free(polynomials);
	if (derivatives != NULL) free(derivatives);
	if (aoverrpowers != NULL) free(aoverrpowers);
	if (times != NULL) free(times);
	if (gnm != NULL) free(gnm);
	if (hnm != NULL) free(hnm);
	if (gnmNow != NULL) free(gnmNow);
	if (hnmNow != NULL) free(hnmNow);
	for (int i = 0; i < NMAGVARS; i++)
	{
		if (magVariables[i] != NULL) free(magVariables[i]);
	}
	if (f != NULL) fclose(f);

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

int calculateField(double r, double theta, double phi, int maxN, double *gnm, double *hnm, double *polynomials, double *derivatives, double *aoverrpowers, double *bn, double *be, double *bc)
{
	double a = 6371.2;
	double aoverr = 1.0;
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
	status = gsl_sf_legendre_deriv_array(GSL_SF_LEGENDRE_SCHMIDT, maxN, cos(theta), polynomials, derivatives);
	if (status)
	{
		printf("status: %s\n", gsl_strerror(status));
		return -1;
	}
	aoverr = a / r;
	aoverrpowers[0] = aoverr * aoverr; // (a/r)^n+1, n starting at 1
	for (int i = 2; i <= maxN; i++)
	{
		aoverrpowers[i-1] = aoverrpowers[i-2] * aoverr;
	}

	for (int n = 1; n <= maxN; n++)
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