/*

    CHAOS: chaos_calc.c

    Copyright (C) 2023  Johnathan K Burchill

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

// Adapted from chaos.c

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

// declare functions
char infoHeader[50] = "";

// Handle Ctrl-C
// From https://stackoverflow.com/questions/4217037/catch-ctrl-c-in-c
sig_atomic_t keep_running = 1;
static void sig_handler(int ignored)
{
    (void)ignored; // Avoid warning
    keep_running = 0;
}

void usage(const char *name);

enum CHAOS_STATUS
{
    CHAOS_STATUS_OK = 0,
    CHAOS_STATUS_POINTERS,
    CHAOS_STATUS_INPUT_FILE,
    CHAOS_STATUS_MEM
};

typedef struct Data
{
    double unixTime;
    double latitude;
    double longitude;
    double altitude;
    double bCoreN;
    double bCoreE;
    double bCoreC;
    double bCrustN;
    double bCrustE;
    double bCrustC;
} Data;

int loadInputsFromFile(char *inFile, Data **data, size_t *nInputs, bool verbose);

int main (int argc, char **argv)
{

	int status = 0;

	// Handle Ctrl-C
	signal(SIGINT, sig_handler);

	time_t processingStartTime = time(NULL);

	ChaosCoefficients coeffs = {0};

	size_t nInputs = 0;
    Data *data = NULL;

    int optionsCount = 0;
    bool overwrite = false;
    bool verbose = false;

	for (int i = 0; i < argc; i++)
	{
		if (strcmp(argv[i], "--about") == 0)
		{
            fprintf(stdout, "chaos_calc - CHAOS core and crustal magnetic field model calculator, version %s.\n", SOFTWARE_VERSION);
            fprintf(stdout, "Copyright (C) 2023  Johnathan K Burchill\n");
            fprintf(stdout, "This program comes with ABSOLUTELY NO WARRANTY.\n");
            fprintf(stdout, "This is free software, and you are welcome to redistribute it\n");
            fprintf(stdout, "under the terms of the GNU General Public License.\n");
			exit(EXIT_SUCCESS);
		}
		else if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--overwrite") == 0)
		{
            optionsCount++;
            overwrite = true;
		}
		else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0)
		{
            optionsCount++;
            verbose = true;
		}
		else if (strcmp(argv[i], "--help") == 0)
		{
			usage(argv[0]);
			exit(EXIT_SUCCESS);
		}
        else if (strncmp(argv[i], "--", 2) == 0)
        {
            fprintf(stderr, "Unknown option %s\n", argv[i]);
            exit(EXIT_FAILURE);
        }
	}

	if ((argc - optionsCount) != 3 || argv[1] == NULL)
	{
		usage(argv[0]);
		exit(EXIT_FAILURE);
	}

	char *inFile = (char*)argv[1];
	char *coeffDir = argv[2];

	char fullOutputFilename[FILENAME_MAX] = {0};
	status = snprintf(fullOutputFilename, FILENAME_MAX-4, "%s.out", inFile);
	if (status < 0)
	{
		fprintf(stderr, "Could not construct full output filename.\n");
		exit(EXIT_FAILURE);
	}
	if (access(fullOutputFilename, F_OK) == 0 && !overwrite)
	{
		printf("Output CDF file exists. Use -f to overwrite. Exiting.\n");
		exit(EXIT_FAILURE);
	}

	status = loadModelCoefficients(coeffDir, &coeffs);
	if (status != SHC_OK || !coeffs.initialized)
	{
		fprintf(stderr, "Error reading CHAOS 7 coefficients files: return code = %d.\n", status);
		goto cleanup;
	}

    status = loadInputsFromFile(inFile, &data, &nInputs, verbose);
    if (status != CHAOS_STATUS_OK)
    {
        fprintf(stderr, "Error loading inputs from %s: return code = %d\n", inFile, status);
        goto cleanup;
    }
    if (nInputs == 0)
    {
        fprintf(stderr, "No inputs found in %s\n", inFile);
        goto cleanup;
    }

    // Date from first entry
    time_t t = (time_t)data[0].unixTime;
    struct tm *date = gmtime(&t);
    if (date == NULL)
    {
        fprintf(stderr, "Error interpreting first input's time.\n");
        goto cleanup;
    }
    int year = date->tm_year + 1900;
    int month = date->tm_mon + 1;
    int day = date->tm_mday;

	// Because MAG inputs are from daily CDF files, we only interpolate once.
	// This is sufficient accuracy for space physics
    status = interpolateSHCCoefficients(&coeffs, year, month, day);
	if (status != SHC_OK)
	{
		fprintf(stderr, "Could not interpolate model coefficients: return code = %d.\n", status);
		goto cleanup;
	}

    // Calculate and print output to file
    Data *p = NULL;
	double degrees = M_PI / 180.0;
	double a = EARTH_RADIUS_KM;
	double inputTime = 0.;
	double r = 0.;
	double theta = 0.0 * degrees;
	double phi = 0.0 * degrees;
    for (int i = 0; i < nInputs && keep_running; i++)
    {
        p = &data[i];
        r = p->altitude + EARTH_RADIUS_KM;
        theta = (90.0 - p->latitude) * degrees;
        phi = p->longitude * degrees;
        status = calculateField(r, theta, phi, &coeffs.core, &p->bCoreN, &p->bCoreE, &p->bCoreC);
        if (status != CHAOS_MODEL_OK)
        {
            fprintf(stderr, "Could not calculate core field: return code = %d\n", status);
            goto cleanup;
        }
        status = calculateField(r, theta, phi, &coeffs.crust, &p->bCrustN, &p->bCrustE, &p->bCrustC);
        if (status != CHAOS_MODEL_OK)
        {
            fprintf(stderr, "Could not calculate crustal field: return code = %d\n", status);
            goto cleanup;
        }
        fprintf(stdout, "%lf %lf %lf %lf %lf %lf %lf\n", p->unixTime, p->latitude, p->longitude, p->altitude, p->bCoreN + p->bCrustN, p->bCoreE + p->bCrustE, p->bCoreC + p->bCrustC);
    }

cleanup:
	freeChaosCoefficients(&coeffs);
    free(data);

	return 0;
}

void usage(const char* name)
{
	printf("Usage: %s <inputfile> <chaosModelCoefficientsDir> [options...]\n", name);
	printf(" <chaosModelCoefficientsDir>: directory containing SHC files\n");
    printf("Options:\n");
    printf(" --overwrite (-f): force overwriting existing .out file if it exists.\n");
    printf(" --verbse (-v): write a little more.\n");
	printf(" --about: print version and license information.\n");
    printf(" --help: print this message.\n");

	return;
}


int loadInputsFromFile(char *inFile, Data **data, size_t *nInputs, bool verbose)
{
    int status = CHAOS_STATUS_OK;

    if (inFile == NULL || data == NULL || nInputs == NULL)
        return CHAOS_STATUS_POINTERS;

    if (verbose)
    	printf("Reading inputs from %s\n", inFile);

    FILE *file = fopen(inFile, "r");
    if (file == NULL)
        return CHAOS_STATUS_INPUT_FILE;

    char buf[256];

    int nLines = 0;
    while (fgets(buf, 256, file) != NULL)
        nLines++;

    *data = calloc(nLines, sizeof **data);
    if (*data == NULL)
    {
        fclose(file);
        return CHAOS_STATUS_MEM;
    }

    printf("nLines: %d\n", nLines);
    status = fseek(file, 0, SEEK_SET);
    if (status != 0)
    {
        fclose(file);
        return CHAOS_STATUS_INPUT_FILE;
    }

    Data *p = NULL;
    int scanned = 0;
    for (int i = 0; i < nLines; i++)
    {
        p = &(*data)[i];
        scanned = fscanf(file, "%lg %lg %lg %lg\n", &p->unixTime, &p->latitude, &p->longitude, &p->altitude);
        if (scanned != 4)
        {
            fclose(file);
            return CHAOS_STATUS_INPUT_FILE;
        }
    }

    fclose(file);

    *nInputs = nLines;

    return CHAOS_STATUS_OK;
}
