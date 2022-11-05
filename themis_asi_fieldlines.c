#include "trace.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <signal.h>
#include <string.h>
#include <unistd.h>

#include <libgen.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <cdf.h>

#include <omp.h>

sig_atomic_t keep_running;
char infoHeader[50];

int main (int argc, char *argv[])
{

    gsl_set_error_handler_off();
    sprintf(infoHeader, "%s", "themis_asi_fieldlines: ");
    keep_running = true;

    int status = CHAOS_TRACE_OK;

    int nOptions = 0;

    double minimumAltitudekm = 0.0;
    double accuracy = 0.001;

    for (int i = 0; i < argc; i++)
    {
        if (strncmp("--minimum-altitude-km=", argv[i], 21) == 0)
        {
            char *lastParsedChar = argv[i]+21;
            double value = strtod(argv[i] + 21, &lastParsedChar);
            if (lastParsedChar == argv[i] + 21)
            {
                fprintf(stderr, "%s: unable to parse %s\n", argv[0], argv[i]);
                exit(EXIT_FAILURE);
            }
            minimumAltitudekm = value;
            nOptions++;
        }
        if (strncmp("--accuracy=", argv[i], 11) == 0)
        {
            char *lastParsedChar = argv[i]+11;
            double value = strtod(argv[i] + 11, &lastParsedChar);
            if (lastParsedChar == argv[i] + 11)
            {
                fprintf(stderr, "%s: unable to parse %s\n", argv[0], argv[i]);
                exit(EXIT_FAILURE);
            }
            accuracy = value;
            nOptions++;
        }
    }


    if (argc - nOptions != 8)
    {
        printf("Incorrect number of arguments.\n");
        printf("usage: %s themisL2File coeffDir year month day startAltkm swarmAltkm [--minimum-altitude-km=value] [--accuracy=value]\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    char *themisL2File = argv[1];
    char *coeffDir = argv[2];
    int year = atoi(argv[3]);
    int month = atoi(argv[4]);
    int day = atoi(argv[5]);
    double startAlt = atof(argv[6]);
    double swarmAlt = atof(argv[7]);

    float finalLatitudes[257][257] = {0.0};
    float finalLongitudes[257][257] = {0.0};
    float finalAltitudeskm[257][257] = {0.0};

    // Open Themis L2 file and get geographic position of each 256x256 pixel for specified altitude. (65535 positions)

    if (access(themisL2File, F_OK) != 0)
    {
        fprintf(stderr, "%s: Unable to access the Themis Level 2 calibration CDF file.\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    CDFid cdf = NULL;
    CDFstatus cdfStatus = CDF_OK;

    cdfStatus = CDFopen(themisL2File, &cdf);
    if (cdfStatus != CDF_OK)
    {
        fprintf(stderr, "%s: Unable to open the Themsi Level 2 calibration CDF file.\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    // Each site name is 4 ASCII characters
    char *base = basename(themisL2File);
    if (base == NULL)
    {
        fprintf(stderr, "%s: Could not get basename of themis level 2 CDF file.\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    char site[5] = {0};
    strncpy(site, base + 11, 4);

    char latVar[255] = {0};
    sprintf(latVar, "thg_asf_%s_glat", site);
    
    char lonVar[255] = {0};
    sprintf(lonVar, "thg_asf_%s_glon", site);

    char altVar[255] = {0};
    sprintf(altVar, "thg_asf_%s_alti", site);


    CDFdata data = NULL;
    long nRecs = 0;
    long dataType = 0;
    long nElem = 0;
    long nDims = 0;
    long dimSizes[CDF_MAX_DIMS] = {0};
    long recVary = 0;
    long dimsVary[CDF_MAX_DIMS] = {0};

    // Altitude
    cdfStatus = CDFreadzVarAllByVarName(cdf, altVar, &nRecs, &dataType, &nElem, &nDims, dimSizes, &recVary, dimsVary, &data);
    if (cdfStatus != CDF_OK)
    {
        fprintf(stderr, "%s: Unable to read altitude data.\n", argv[0]);
        CDFclose(cdf);
        exit(EXIT_FAILURE);
    }

    // 3 floats
    float altitudes[3] = {0.0};
    bool validAltitude = false;
    int altitudeIndex = -1;
    for (int i = 0; i < 3; i++)
    {
        altitudes[i] = ((float*)data)[i] / 1000.0;
        if ((int)round(altitudes[i]) == (int)round(startAlt))
        {
            validAltitude = true;
            altitudeIndex = i;
        }
    }
    CDFdataFree(data);

    if (validAltitude != true)
    {
        fprintf(stderr, "invalid startAlt. Use one of (in km):\n");
        for (int i = 0; i < 3; i++)
        {
            fprintf(stderr, "\t%.0f\n", altitudes[i]);
        }
        CDFclose(cdf);
        exit(EXIT_FAILURE);
    }

    // 257 x 257 latitudes and longitudes at corners of pixels
    float *latitudes = (float*)calloc(257 * 257, sizeof(float));
    float *longitudes = (float*)calloc(257 * 257,  sizeof(float));
    if (latitudes == NULL || longitudes == NULL)
    {
        fprintf(stderr, "%s: Could not allocate memory.\n", argv[0]);
        CDFclose(cdf);
        exit(EXIT_FAILURE);        
    }

    // Epochs
    cdfStatus = CDFreadzVarAllByVarName(cdf, "range_epoch", &nRecs, &dataType, &nElem, &nDims, dimSizes, &recVary, dimsVary, &data);

    // Latitudes at corners of pixels
    cdfStatus = CDFreadzVarAllByVarName(cdf, latVar, &nRecs, &dataType, &nElem, &nDims, dimSizes, &recVary, dimsVary, &data);
    if (cdfStatus != CDF_OK)
    {
        fprintf(stderr, "%s: Unable to read latitude data.\n", argv[0]);
        CDFclose(cdf);
        exit(EXIT_FAILURE);
    }
    for (int i = altitudeIndex; i < 257*257*3; i+=3)
    {
        latitudes[i/3] = ((float*)data)[(257*257*3) + i];
        // latitudes[i/3] = ((float*)data)[i];
    }
    CDFdataFree(data);

    cdfStatus = CDFreadzVarAllByVarName(cdf, lonVar, &nRecs, &dataType, &nElem, &nDims, dimSizes, &recVary, dimsVary, &data);
    if (cdfStatus != CDF_OK)
    {
        fprintf(stderr, "%s: Unable to read longitude data.\n", argv[0]);
        CDFclose(cdf);
        exit(EXIT_FAILURE);
    }
    for (int i = altitudeIndex; i < 257*257*3; i+=3)
    {
        longitudes[i/3] = ((float*)data)[(257*257*3) + i];
        // longitudes[i/3] = ((float*)data)[i];
    }
    CDFdataFree(data);

    // for (int i = 0; i < 257*257; i++)
    // {
    //     printf("%.1f %.1f\n", latitudes[i], longitudes[i]);
    // }

	ChaosCoefficients coeffs = {0};
    status = initializeTracer(coeffDir, year, month, day, &coeffs);
    long steps = 0;

    double latitude = 0.0;
    double longitude = 0.0;
    double altitude = 0.0;
    double lat = 0.0;
    double lon = 0.0;



    for (int i = 0; i < 257; i++)
    {
        if (i % 26 == 0)
            fprintf(stderr, "i = %d of 257\n", i);
#pragma omp parallel for private(lat,lon,status,latitude,longitude,altitude,steps)
        for (int j = 0; j < 257; j++)
        {
            lat = latitudes[i*257+j];
            lon = longitudes[i*257+j];
            if (!isfinite(lat) || !isfinite(lon))
                continue;
            status = trace(&coeffs, -1, accuracy, lat, lon, startAlt, minimumAltitudekm, swarmAlt, &latitude, &longitude, &altitude, &steps);
            finalLatitudes[i][j] = latitude;
            finalLongitudes[i][j] = longitude;
            finalAltitudeskm[i][j] = altitude;
        }


    }

    for (int i = 0; i < 257; i++)
    {
        for (int j = 0; j < 257; j++)
        {
            lat = latitudes[i*257+j];
            lon = longitudes[i*257+j];
            if (!isfinite(lat) || !isfinite(lon))
                continue;
            printf("%d %d %lf %lf %lf\n", i, j, finalLatitudes[i][j], finalLongitudes[i][j], finalAltitudeskm[i][j]*1000.0);
        }


    }

    CDFclose(cdf);
    free(latitudes);
    free(longitudes);
    freeChaosCoefficients(&coeffs);

    return EXIT_SUCCESS;
}

