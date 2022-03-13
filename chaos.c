#include "cdf_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_math.h>

#define NMAGVARS 5

void usage(const char *name);

char infoHeader[50];

int main (int argc, char **argv)
{

	int status = 0;

	printf("Hi\n");

	double degrees = M_PI / 180.0;
	double a = 6371.2;
	double time = 0.;
	double r = 0.;
	double theta = 0.0 * degrees;
	double phi = 0.0 * degrees;

	double *polynomials = NULL;
	double *aoverrpowers = NULL;
	double *times = NULL;
	double *gnm = NULL;
	double *hnm = NULL;

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
	aoverrpowers = malloc(maxN * sizeof(double));
	times = (double*)malloc(nTimes * sizeof(double));
	// Uses more memory than needed. Could be revised
	gnm = (double*)malloc(maxN * (maxN + 2) * nTimes * sizeof(double));
	hnm = (double*)malloc(maxN * (maxN + 2) * nTimes * sizeof(double));
	if (polynomials == NULL || aoverrpowers == NULL || times == NULL || gnm == NULL || hnm == NULL)
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

	char *magVariableNames[NMAGVARS] = {
		"Timestamp",
		"Latitude",
		"Longitude",
		"Radius",
		"B_NEC"
	};
	// Magnetic field input data

	if (getInputFilename(satellite, year, month, day, magDir, "HR_1B", magFilename))
    {
        fprintf(stdout, "%sMAG LR_1B input file is not available. Exiting.\n", infoHeader);
        exit(1);
    }
	loadCdf(magFilename, magVariableNames, NMAGVARS, magVariables, &nInputs);
	if (nInputs == 0)
	{
		printf("Found no inputs from MAG_HR file.\n");
		goto cleanup;
	}

	double x = 0;
	double magneticPotential = 0.0;
	double magneticPotentialN = 0.0;

	double aoverr = 1.0;

	printf("Calculating magnetic potential...\n");
	for (size_t t = 0; t < nInputs; t+=25)
	{
		// Interpolate coefficients to magnetic time.
		time = ((double*)magVariables[0])[t];
	 	theta = (90.0 - ((double*)magVariables[1])[t]) * degrees;
		phi = ((double*)magVariables[2])[t] * degrees;
		r = ((double*)magVariables[3])[t]/1000.;
		status = gsl_sf_legendre_array(GSL_SF_LEGENDRE_SCHMIDT, maxN, cos(theta), polynomials);
		if (status)
		{
			printf("status: %s\n", gsl_strerror(status));
			goto cleanup;
		}
		aoverr = a / r;
		aoverrpowers[0] = aoverr * aoverr; // (a/r)^n+1, n starting at 1
		for (int i = 2; i <= maxN; i++)
		{
			aoverrpowers[i-1] = aoverrpowers[i-2] * aoverr;
		}

		magneticPotential = 0.0;
		gRead = 0;
		hRead = 0;
		for (int n = 1; n <= maxN; n++)
		{
			magneticPotentialN = gnm[gRead*nTimes] * polynomials[gsl_sf_legendre_array_index(n, 0)];
			gRead++;
			for (int m = 1; m <= n; m++)
			{
				magneticPotentialN += (gnm[gRead*nTimes + 0] * cos(phi) + hnm[hRead*nTimes + 0] * sin(phi) ) * polynomials[gsl_sf_legendre_array_index(n, m)];
				gRead++;
				hRead++;
			}
			magneticPotential = magneticPotentialN * aoverrpowers[n-1] * a;
		}
		printf("magneticPotential(time=%lf, r=%lf, colatitude=%lf, longitude=%lf) = %lf\n", time, r, theta/degrees, phi/degrees, magneticPotential);
	}

cleanup:
	if (polynomials != NULL) free(polynomials);
	if (aoverrpowers != NULL) free(aoverrpowers);
	if (times != NULL) free(times);
	if (gnm != NULL) free(gnm);
	if (hnm != NULL) free(hnm);
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