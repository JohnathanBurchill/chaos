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
#include "shc.h"
#include "model.h"
#include "chaos_settings.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <signal.h>

// #include <gsl/gsl_errno.h>
// #include <gsl/gsl_sf_legendre.h>
// #include <gsl/gsl_math.h>

#define NMAGVARS 5

// declare functions

// Handle Ctrl-C
// From https://stackoverflow.com/questions/4217037/catch-ctrl-c-in-c
sig_atomic_t keep_running = 1;
static void sig_handler(int ignored)
{
    (void)ignored; // Avoid warning
    keep_running = 0;
}

void usage(const char *name);

char infoHeader[50];

int main (int argc, char **argv)
{

	int status = 0;

	// Handle Ctrl-C
	signal(SIGINT, sig_handler);

	time_t processingStartTime = time(NULL);

	char magFilename[FILENAME_MAX];
	char outputFilename[FILENAME_MAX];
	double beginTime = 0.0;
	double endTime = 0.0;

	ChaosCoefficients coeffs = {0};

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

	for (int i = 0; i < argc; i++)
	{
		if (strcmp(argv[i], "--about") == 0)
		{
			// Adapted from TII Ion Drift Processor
            fprintf(stdout, "chaos - CHAOS core and crustal magnetic field model and residual field processor, version %s.\n", SOFTWARE_VERSION);
            fprintf(stdout, "Copyright (C) 2022  Johnathan K Burchill\n");
            fprintf(stdout, "This program comes with ABSOLUTELY NO WARRANTY.\n");
            fprintf(stdout, "This is free software, and you are welcome to redistribute it\n");
            fprintf(stdout, "under the terms of the GNU General Public License.\n");
			exit(EXIT_SUCCESS);
		}
		else if (strcmp(argv[i], "--help") == 0)
		{
			usage(argv[0]);
			exit(EXIT_SUCCESS);
		}
	}

	if (argc != 6)
	{
		usage(argv[0]);
		exit(EXIT_FAILURE);
	}

	char *satDate = argv[1];
	char *magDataset = argv[2];
	char *coeffDir = argv[3];
	char *magDir = argv[4];
	char *outputDir = argv[5];

	if (strcmp(magDataset, "LR_1B") != 0 && strcmp(magDataset, "HR_1B") != 0)
	{
		fprintf(stderr, "Expected 'LR_1B' or 'HR_1B' for magDataset.\n");
		usage(argv[0]);
		exit(EXIT_FAILURE);
	}

	char satellite = satDate[0];
	long year, month, day;
    int valuesRead = sscanf(satDate+1, "%4ld%2ld%2ld", &year, &month, &day);
    if (valuesRead != 3)
    {
        fprintf(stderr, "CHAOS called as:\n \"%s %s %s %s %s %s\"\n Unable to parse date \"\". Exiting.\n", argv[0], argv[1], argv[2], argv[3], argv[4], argv[5]);
		exit(EXIT_FAILURE);
    }
	sprintf(infoHeader, "CHAOS %s: ", satDate);

	status = getOutputFilename(satellite, year, month, day, outputDir, &beginTime, &endTime, outputFilename, magDataset);
	if (status != 0)
	{
		fprintf(stderr, "Could not get output filename.\n");
		exit(EXIT_FAILURE);
	}
	char fullOutputFilename[FILENAME_MAX] = {0};
	status = snprintf(fullOutputFilename, FILENAME_MAX-4, "%s.cdf", outputFilename);
	if (status < 0)
	{
		fprintf(stderr, "Could not construct full output filename.\n");
		exit(EXIT_FAILURE);
	}
	if (access(fullOutputFilename, F_OK) == 0)
	{
		printf("%sOutput CDF file exists. Exiting.\n", infoHeader);
		exit(EXIT_FAILURE);
	}

	// Interpolate samples every interpolationSkip
	// 4 s
	int interpolationSkip = 4;
	if (strcmp(magDataset, "HR_1B") == 0)
		interpolationSkip = 200;


	status = loadModelCoefficients(coeffDir, &coeffs);
	if (status != SHC_OK || !coeffs.initialized)
	{
		fprintf(stderr, "%sError reading CHAOS 7 coefficients files: return code = %d.\n", infoHeader, status);
		goto cleanup;
	}

	// Because MAG inputs are from daily CDF files, we only interpolate once.
	// This is sufficient accuracy for space physics
	status = interpolateSHCCoefficients(&coeffs, year, month, day);
	if (status != SHC_OK)
	{
		fprintf(stderr, "%sCould not interpolate model coefficients: return code = %d.\n", infoHeader, status);
		goto cleanup;
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
	if (getInputFilename(satellite, year, month, day, magDir, magDataset, magFilename))
	// if (getInputFilename(satellite, year, month, day, magDir, "HR_1B", magFilename))
    {
        fprintf(stdout, "%sMAG input file is not available. Exiting.\n", infoHeader);
        exit(1);
    }
	printf("%sReading inputs from %s\n", infoHeader, magFilename);
	loadCdf(magFilename, magVariableNames, NMAGVARS, magVariables, &nInputs);
	if (nInputs == 0)
	{
		fprintf(stderr, "%sFound no measurements in MAG file.\n", infoHeader);
		goto cleanup;
	}

	// Measured fields
	dbMeas = (double*)malloc(nInputs * 3 * sizeof(double));

	// Model fields
	bCore = (double*)malloc(nInputs * 3 * sizeof(double));
	bCrust = (double*)malloc(nInputs * 3 * sizeof(double));
	if (dbMeas == NULL || bCore == NULL || bCrust == NULL)
	{
		fprintf(stderr, "%sMemory issue.\n", infoHeader);
		goto cleanup;
	}

	status = calculateResiduals(&coeffs, interpolationSkip, magVariables, nInputs, bCore, bCrust, dbMeas);
	if (status != CHAOS_MODEL_OK)
	{
		fprintf(stderr, "%sCould not calculate all residuals: return code = %d\n", infoHeader, status);
		goto cleanup;
	}

	if (keep_running == 0)
	{
		fprintf(stderr, "%sInterrupted (SIGINT).\n", infoHeader);
		goto cleanup;
	}

	status = exportCdf(outputFilename, magFilename, &coeffs, satellite, magDataset, EXPORT_VERSION_STRING, (double*)magVariables[0], (double*)magVariables[1], (double*)magVariables[2], (double*)magVariables[3], bCore, bCrust, dbMeas, nInputs);
	if (status != 0)
	{
		fprintf(stderr, "%sCould not export fields: return code = %d\n", infoHeader, status);
	}

	// time_t processingStopTime = time(NULL);
	// exportMetaInfo(outputFilename, magFilename, coreFile, crustalFile, nInputs, processingStartTime, processingStopTime);

cleanup:
	freeChaosCoefficients(&coeffs);

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
	printf("Usage: %s XYYYYMMDD magDataset chaosModelCoefficientsDir magCdfDir outputDir\n", name);
	printf(" X: satellite letter A, B, or C\n");
	printf(" YYYYMMDD: year, month, day\n");
	printf(" magDataset:\n");
	printf("\tLR_1B: Input files are  1 Hz data.\n");
	printf("\tHR_1B: Input files are 50 Hz data.\n");
	printf(" chaosModelCoefficientsDir: directory containing SHC files\n");
	printf(" magCdfDir: directory containing MAG_HR CDFs\n");
	printf(" outputDir: directory to store magnetic field vectors\n");

	return;
}


