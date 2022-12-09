#include "trace.h"
#include "model.h"
#include "util.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <signal.h>
#include <string.h>
#include <libgen.h>
#include <unistd.h>
#include <time.h>

#include <libgen.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <cdf.h>

sig_atomic_t keep_running;
char infoHeader[50];

#define PROGRAM_VERSION_STRING "20221208"
#define IMAGE_COLUMNS 256
#define IMAGE_ROWS 256

int addVariableAttributes(CDFid cdf, char *name, char *description, char *units);

int main (int argc, char *argv[])
{

    gsl_set_error_handler_off();
    sprintf(infoHeader, "%s", "themis_asi_fieldlines: ");
    keep_running = true;

    double processingStart = 0.0;
    double unixTime = (double) time(NULL);
    UnixTimetoEPOCH(&unixTime, &processingStart, 1);


    int status = CHAOS_TRACE_OK;

    int nOptions = 0;

    double minimumAltitudekm = 0.0;
    double accuracy = 0.001;

    bool showProgress = false;

    char *exportDir = ".";

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
        if (strcmp("--progress", argv[i]) == 0)
        {
            nOptions++;
            showProgress = true;
        }
        if (strncmp("--exportdir=", argv[i], 12) == 0)
        {
            nOptions++;
            exportDir = argv[i] + 12;
        }
    }

    if (argc - nOptions != 5)
    {
        printf("Incorrect number of arguments.\n");
        printf("usage: %s calibrationFile coeffDir startAltkm targetAltkm [--minimum-altitude-km=value] [--accuracy=value] [--progress]\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    char *calibrationFile = argv[1];
    char *coeffDir = argv[2];
    double startAltKm = atof(argv[3]);
    double targetAltKm = atof(argv[4]);

    // Open Themis L2 file and get geographic position of each 256x256 pixel for specified altitude. (65535 positions)

    if (access(calibrationFile, F_OK) != 0)
    {
        fprintf(stderr, "%s: Unable to access the calibration CDF file.\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    CDFid cdf = NULL;
    CDFstatus cdfStatus = CDF_OK;

    cdfStatus = CDFopen(calibrationFile, &cdf);
    if (cdfStatus != CDF_OK)
    {
        fprintf(stderr, "%s: Unable to open the calibration CDF file.\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    char siteAbbreviation[CDF_ATTR_NAME_LEN256+1] = {0};
    cdfStatus = CDFgetzVarRecordData(cdf, CDFgetVarNum(cdf, "SiteAbbreviation"), 0, siteAbbreviation);
    if (cdfStatus != CDF_OK)
    {
        fprintf(stderr, "%s: Unable to read site abbreviation from CDF file.\n", argv[0]);
        CDFclose(cdf);
        exit(EXIT_FAILURE);
    }

    if (strlen(siteAbbreviation) != 4)
    {
        fprintf(stderr, "%s: Expected 4-letter site abbreviation, got %s.", argv[0], siteAbbreviation);
        CDFclose(cdf);
        exit(EXIT_FAILURE);
    }

    CDFdata data = NULL;
    long nRecs = 0;
    long dataType = 0;
    long nElem = 0;
    long nDims = 0;
    long dimSizes[CDF_MAX_DIMS] = {0};
    long recVary = 0;
    long dimsVary[CDF_MAX_DIMS] = {0};

    // Site Altitude (geodetic)

    float siteLatDeg = 0.0;
    float siteLonDeg = 0.0;
    float siteAltM = 0.0;

    float elevations[IMAGE_COLUMNS][IMAGE_ROWS];
    float azimuths[IMAGE_COLUMNS][IMAGE_ROWS];


    cdfStatus = CDFgetzVarRecordData(cdf, CDFgetVarNum(cdf, "SiteLatitude"), 0, &siteLatDeg);
    if (cdfStatus != CDF_OK)
    {
        fprintf(stderr, "%s: Unable to read site latitude from CDF file.\n", argv[0]);
        CDFclose(cdf);
        exit(EXIT_FAILURE);
    }
    cdfStatus = CDFgetzVarRecordData(cdf, CDFgetVarNum(cdf, "SiteLongitude"), 0, &siteLonDeg);
    if (cdfStatus != CDF_OK)
    {
        fprintf(stderr, "%s: Unable to read site longitude from CDF file.\n", argv[0]);
        CDFclose(cdf);
        exit(EXIT_FAILURE);
    }
    cdfStatus = CDFgetzVarRecordData(cdf, CDFgetVarNum(cdf, "SiteAltitude"), 0, &siteAltM);
    if (cdfStatus != CDF_OK)
    {
        fprintf(stderr, "%s: Unable to read site altitude from CDF file.\n", argv[0]);
        CDFclose(cdf);
        exit(EXIT_FAILURE);
    }

    double calibrationEpoch = 0.0;
    cdfStatus = CDFgetzVarRecordData(cdf, CDFgetVarNum(cdf, "CalibrationEpoch"), 0, &calibrationEpoch);
    if (cdfStatus != CDF_OK)
    {
        fprintf(stderr, "%s: Unable to read calibration epoch from CDF file.\n", argv[0]);
        CDFclose(cdf);
        exit(EXIT_FAILURE);
    }

    cdfStatus = CDFreadzVarAllByVarName(cdf, "CalibratedElevations", &nRecs, &dataType, &nElem, &nDims, dimSizes, &recVary, dimsVary, &data);
    if (cdfStatus != CDF_OK)
    {
        fprintf(stderr, "%s: Unable to read elevation data from CDF.\n", argv[0]);
        CDFclose(cdf);
        exit(EXIT_FAILURE);
    }
    memcpy(elevations, data, sizeof(float) * IMAGE_COLUMNS * IMAGE_ROWS);
    CDFdataFree(data);

    cdfStatus = CDFreadzVarAllByVarName(cdf, "CalibratedAzimuths", &nRecs, &dataType, &nElem, &nDims, dimSizes, &recVary, dimsVary, &data);
    if (cdfStatus != CDF_OK)
    {
        fprintf(stderr, "%s: Unable to read azimuth data from CDF.\n", argv[0]);
        CDFclose(cdf);
        exit(EXIT_FAILURE);
    }
    memcpy(azimuths, data, sizeof(float) * IMAGE_COLUMNS * IMAGE_ROWS);
    CDFdataFree(data);

    CDFclose(cdf);

    // Check for existing export CDF before getting into the calculations
    long year = 0;
    long month = 0;
    long day = 0;
    long hour = 0;
    long minute = 0;
    long second = 0;
    long millisecond = 0;
    EPOCHbreakdown(calibrationEpoch, &year, &month, &day, &hour, &minute, &second, &millisecond);

    char traceCdf[CDF_PATHNAME_LEN + 1] = {0};
    snprintf(traceCdf, CDF_PATHNAME_LEN + 1, "%s/themis_%s_chaos7_trace_%4ld%02ld%02ldT%02ld%02ld%02ld_%.0lfkm_to_%.0lfkm", exportDir, siteAbbreviation, year, month, day, hour, minute, second, startAltKm, targetAltKm);

    char traceFullCdf[CDF_PATHNAME_LEN + 5] = {0};
    snprintf(traceFullCdf, CDF_PATHNAME_LEN + 5, "%s.cdf", traceCdf);

    if (access(traceFullCdf, F_OK) == 0)
    {
        fprintf(stderr, "Will not overwrite existing CDF. Delete %s and try again.\n", traceFullCdf);
        return EXIT_FAILURE;
    }

    // Using double here for magnetic field line tracing
    double enuCorners[IMAGE_COLUMNS+1][IMAGE_ROWS+1][3] = {0};
    // Geocentric lat, lon, radius
    double geocentricPositionCorners[IMAGE_COLUMNS + 1][IMAGE_ROWS + 1][3];
    double tracedGeocentricPositionCorners[IMAGE_COLUMNS + 1][IMAGE_ROWS + 1][3];

    // Init to NAN. No need to calculate corners at edge of image, where 
    // elevation and azimuth values are NAN anyway.
    for (int c = 0; c < IMAGE_COLUMNS + 1; c++)
        for (int r = 0; r < IMAGE_ROWS + 1; r++)
            for (int k = 0; k < 3; k++)
            {
                enuCorners[c][r][k] = nan("");
                geocentricPositionCorners[c][r][k] = nan("");
                tracedGeocentricPositionCorners[c][r][k] = nan("");
            }

    double enu1[3] = {0};
    double enu2[3] = {0};
    double enu3[3] = {0};
    double enu4[3] = {0};

    double sitePositionXyz[3] = {0};
    geodeticToGeocentric(siteLatDeg, siteLonDeg, siteAltM, NULL, NULL, sitePositionXyz);

    // Values at edge of image are NaNs, no need to handle corners from the two image sides
    for (int c = 1; c < IMAGE_COLUMNS; c++)
    {
        for (int r = 1; r < IMAGE_ROWS; r++)
        {
            elazToEnu((double)elevations[c-1][r-1], (double)azimuths[c-1][r-1], enu1);
            elazToEnu((double)elevations[c][r-1], (double)azimuths[c][r-1], enu2);
            elazToEnu((double)elevations[c-1][r], (double)azimuths[c-1][r], enu3);
            elazToEnu((double)elevations[c][r], (double)azimuths[c][r], enu4);

            if (!isfinite(enu1[0]) || !isfinite(enu2[0]) || !isfinite(enu3[0]) || !isfinite(enu4[0]))
                continue;

            for (int k = 0; k < 3; k++)
                enuCorners[c][r][k] = (double)((enu1[k] + enu2[k] + enu3[k] + enu4[k]) / 4.0);
                            
            // printf("enuCorners[%d][%d] = (%f, %f, %f)\n", c, r, enuCorners[c][r][0], enuCorners[c][r][1], enuCorners[c][r][2]);
    
            // Geocentric latitudes and longitudes at corners of pixels
            lookDirectionToPosition(sitePositionXyz, enuCorners[c][r], startAltKm, geocentricPositionCorners[c][r]);
            // printf("%lf %lf %lf\n", geocentricPositionCorners[c][r][0], geocentricPositionCorners[c][r][1], geocentricPositionCorners[c][r][2]);

        }

    }

	ChaosCoefficients coeffs = {0};


    status = initializeTracer(coeffDir, (int)year, (int)month, (int)day, &coeffs);
    long steps = 0;

    double latitude = 0.0;
    double longitude = 0.0;
    double altitude = 0.0;
    double sphericalAltKm = 0.0;

    printf("MinAltkm: %lf\n", minimumAltitudekm);
    printf("targetAltKm: %lf\n", targetAltKm);
    int progressCounter = 0;
    int expectedTraces = (IMAGE_COLUMNS + 1) * (IMAGE_ROWS + 1);
    for (int i = 0; i < IMAGE_COLUMNS + 1; i++)
    {
        for (int j = 0; j < IMAGE_ROWS + 1; j++)
        {
            progressCounter++;
            if (progressCounter >= 350)
                break;
            if (!isfinite(geocentricPositionCorners[i][j][0]))
                continue;
            sphericalAltKm = geocentricPositionCorners[i][j][2] / 1000.0 - EARTH_RADIUS_KM;
            // Results in geocentric latitude, longitude, and spherical altitude in km(geocentric radius minus mean earth radius)
            status = trace(&coeffs, -1, accuracy, geocentricPositionCorners[i][j][0], geocentricPositionCorners[i][j][1], sphericalAltKm, minimumAltitudekm, targetAltKm, &latitude, &longitude, &altitude, &steps);
            tracedGeocentricPositionCorners[i][j][0] = latitude;
            tracedGeocentricPositionCorners[i][j][1] = longitude;
            tracedGeocentricPositionCorners[i][j][2] = 1000.0 * (altitude + EARTH_RADIUS_KM);
            // printf("pixelCornerIndex(%d, %d): %.2lf -> %.2lf N %.2lf -> %.2lf E %.1lf -> %.1lf km\n", i, j, geocentricPositionCorners[i][j][0], latitude, geocentricPositionCorners[i][j][1], longitude, sphericalAltKm, altitude);
            if (showProgress)
                fprintf(stderr, "\r%6d of %d field lines traced", progressCounter, expectedTraces);

        }

    }

    if (showProgress)
        fprintf(stderr, "\n");

    freeChaosCoefficients(&coeffs);

    // Export trace results to CDF
    cdf = NULL;
    cdfStatus = CDFcreateCDF(traceCdf, &cdf);
    if (cdfStatus != CDF_OK)
    {
        fprintf(stderr, "Could not create CDF\n");
        return EXIT_FAILURE;
    }

    cdfStatus = CDFsetEncoding(cdf, NETWORK_ENCODING);
    if (cdfStatus != CDF_OK)
    {
        fprintf(stderr, "Could not set CDF encoding.\n");
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }

    nDims = 2;
    recVary = VARY;
    dimSizes[0] = IMAGE_COLUMNS + 1;
    dimSizes[1] = IMAGE_ROWS + 1;
    dimsVary[0] = VARY;
    dimsVary[1] = VARY;
    long varNum = 0;
    double outputBuffer[(IMAGE_COLUMNS + 1) * (IMAGE_ROWS + 1)] = {0};

    nDims = 0;
    cdfStatus = CDFcreatezVar(cdf, "SiteAbbreviation", CDF_CHAR, strlen(siteAbbreviation), nDims, dimSizes, recVary, dimsVary, &varNum);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    cdfStatus = CDFputzVarAllRecordsByVarID(cdf, varNum, 1, siteAbbreviation);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }

    nDims = 0;
    cdfStatus = CDFcreatezVar(cdf, "SiteLatitude", CDF_REAL4, 1, nDims, dimSizes, recVary, dimsVary, &varNum);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    cdfStatus = CDFputzVarAllRecordsByVarID(cdf, varNum, 1, &siteLatDeg);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }

    nDims = 0;
    cdfStatus = CDFcreatezVar(cdf, "SiteLongitude", CDF_REAL4, 1, nDims, dimSizes, recVary, dimsVary, &varNum);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    cdfStatus = CDFputzVarAllRecordsByVarID(cdf, varNum, 1, &siteLonDeg);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }

    nDims = 0;
    cdfStatus = CDFcreatezVar(cdf, "SiteAltitude", CDF_REAL4, 1, nDims, dimSizes, recVary, dimsVary, &varNum);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    cdfStatus = CDFputzVarAllRecordsByVarID(cdf, varNum, 1, &siteAltM);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }

    nDims = 2;
    cdfStatus = CDFcreatezVar(cdf, "EmissionLatitudes", CDF_REAL8, 1, nDims, dimSizes, recVary, dimsVary, &varNum);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    int counter = 0;
    for (int i = 0; i < IMAGE_COLUMNS + 1; i++)
        for (int j = 0; j < IMAGE_ROWS + 1; j++)
            outputBuffer[counter++] = geocentricPositionCorners[i][j][0];
    cdfStatus = CDFputzVarAllRecordsByVarID(cdf, varNum, 1, outputBuffer);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }

    cdfStatus = CDFcreatezVar(cdf, "EmissionLongitudes", CDF_REAL8, 1, nDims, dimSizes, recVary, dimsVary, &varNum);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    counter = 0;
    for (int i = 0; i < IMAGE_COLUMNS + 1; i++)
        for (int j = 0; j < IMAGE_ROWS + 1; j++)
            outputBuffer[counter++] = geocentricPositionCorners[i][j][1];
    cdfStatus = CDFputzVarAllRecordsByVarID(cdf, varNum, 1, outputBuffer);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }

    cdfStatus = CDFcreatezVar(cdf, "EmissionRadii", CDF_REAL8, 1, nDims, dimSizes, recVary, dimsVary, &varNum);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    counter = 0;
    for (int i = 0; i < IMAGE_COLUMNS + 1; i++)
        for (int j = 0; j < IMAGE_ROWS + 1; j++)
            outputBuffer[counter++] = geocentricPositionCorners[i][j][2];
    cdfStatus = CDFputzVarAllRecordsByVarID(cdf, varNum, 1, outputBuffer);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }

    nDims = 0;
    cdfStatus = CDFcreatezVar(cdf, "EmissionAltitudeGeodetic", CDF_REAL8, 1, nDims, dimSizes, recVary, dimsVary, &varNum);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    double emissionAltitude = 1000.0 * startAltKm;
    cdfStatus = CDFputzVarAllRecordsByVarID(cdf, varNum, 1, &emissionAltitude);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }

    nDims = 2;
    cdfStatus = CDFcreatezVar(cdf, "TargetLatitudes", CDF_REAL8, 1, nDims, dimSizes, recVary, dimsVary, &varNum);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    counter = 0; 
    for (int i = 0; i < IMAGE_COLUMNS + 1; i++)
        for (int j = 0; j < IMAGE_ROWS + 1; j++)
            outputBuffer[counter++] = tracedGeocentricPositionCorners[i][j][0];
    cdfStatus = CDFputzVarAllRecordsByVarID(cdf, varNum, 1, outputBuffer);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }

    cdfStatus = CDFcreatezVar(cdf, "TargetLongitudes", CDF_REAL8, 1, nDims, dimSizes, recVary, dimsVary, &varNum);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    counter = 0;
    for (int i = 0; i < IMAGE_COLUMNS + 1; i++)
        for (int j = 0; j < IMAGE_ROWS + 1; j++)
            outputBuffer[counter++] = tracedGeocentricPositionCorners[i][j][1];
    cdfStatus = CDFputzVarAllRecordsByVarID(cdf, varNum, 1, outputBuffer);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }

    cdfStatus = CDFcreatezVar(cdf, "TargetRadii", CDF_REAL8, 1, nDims, dimSizes, recVary, dimsVary, &varNum);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    counter = 0;
    for (int i = 0; i < IMAGE_COLUMNS + 1; i++)
        for (int j = 0; j < IMAGE_ROWS + 1; j++)
            outputBuffer[counter++] = tracedGeocentricPositionCorners[i][j][2];
    cdfStatus = CDFputzVarAllRecordsByVarID(cdf, varNum, 1, outputBuffer);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }

    nDims = 0;
    cdfStatus = CDFcreatezVar(cdf, "TargetAltitudeGeodetic", CDF_REAL8, 1, nDims, dimSizes, recVary, dimsVary, &varNum);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    double targetAltitude = 1000.0 * targetAltKm;
    cdfStatus = CDFputzVarAllRecordsByVarID(cdf, varNum, 1, &targetAltitude);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }

    nDims = 0;
    cdfStatus = CDFcreatezVar(cdf, "CalibrationEpoch", CDF_EPOCH, 1, nDims, dimSizes, recVary, dimsVary, &varNum);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    cdfStatus = CDFputzVarAllRecordsByVarID(cdf, varNum, 1, &calibrationEpoch);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }

    // Global attributes
    long attrNum = 0;
    long entry = 0;
    cdfStatus = CDFcreateAttr(cdf, "Author", GLOBAL_SCOPE, &attrNum);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    char *name = "Johnathan K. Burchill";
    cdfStatus = CDFputAttrgEntry(cdf, attrNum, entry, CDF_CHAR, strlen(name), name);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }

    cdfStatus = CDFcreateAttr(cdf, "Mission", GLOBAL_SCOPE, &attrNum);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    char *mission = "THEMIS ASI";
    cdfStatus = CDFputAttrgEntry(cdf, attrNum, entry, CDF_CHAR, strlen(mission), mission);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }

    cdfStatus = CDFcreateAttr(cdf, "Site", GLOBAL_SCOPE, &attrNum);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    cdfStatus = CDFputAttrgEntry(cdf, attrNum, entry, CDF_CHAR, strlen(siteAbbreviation), siteAbbreviation);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }

    cdfStatus = CDFcreateAttr(cdf, "Processor", GLOBAL_SCOPE, &attrNum);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    char *processorString = "themis_asi_fieldlines version " PROGRAM_VERSION_STRING;
    cdfStatus = CDFputAttrgEntry(cdf, attrNum, entry, CDF_CHAR, strlen(processorString), processorString);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }

    cdfStatus = CDFcreateAttr(cdf, "ProcessingStart", GLOBAL_SCOPE, &attrNum);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    char date[EPOCHx_STRING_MAX+1];
    char attrformat[EPOCHx_FORMAT_MAX+1] = "UTC=<year>-<mm.02>-<dom.02>T<hour>:<min>:<sec>";
    encodeEPOCHx(processingStart, attrformat, date);
    cdfStatus = CDFputAttrgEntry(cdf, attrNum, entry, CDF_CHAR, strlen(date), date);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }

    cdfStatus = CDFcreateAttr(cdf, "ProcessingStop", GLOBAL_SCOPE, &attrNum);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    double processingStop = 0.0;
    unixTime = (double) time(NULL);
    UnixTimetoEPOCH(&unixTime, &processingStop, 1);
    encodeEPOCHx(processingStop, attrformat, date);
    cdfStatus = CDFputAttrgEntry(cdf, attrNum, entry, CDF_CHAR, strlen(date), date);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }

    cdfStatus = CDFcreateAttr(cdf, "ProcessingCommand", GLOBAL_SCOPE, &attrNum);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    size_t commandLength = 0;
    for (int i = 0; i < argc; i++)
        commandLength += strlen(argv[i]) + 1;
    commandLength--;
    char *command = calloc(commandLength, 1);
    if (command == NULL)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    for (int i = 0; i < argc; i++)
    {
        command = strncat(command, argv[i], commandLength);
        if (i != argc-1)
            command = strcat(command, " ");
    }
    cdfStatus = CDFputAttrgEntry(cdf, attrNum, entry, CDF_CHAR, strlen(command), command);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    free(command);

    cdfStatus = CDFcreateAttr(cdf, "Filename", GLOBAL_SCOPE, &attrNum);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    char *cdfBasename = basename(traceFullCdf);
    cdfStatus = CDFputAttrgEntry(cdf, attrNum, entry, CDF_CHAR, strlen(cdfBasename), cdfBasename);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }

    cdfStatus = CDFcreateAttr(cdf, "ReferenceCalibrationFilename", GLOBAL_SCOPE, &attrNum);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    char *calFile = basename(calibrationFile);
    cdfStatus = CDFputAttrgEntry(cdf, attrNum, entry, CDF_CHAR, strlen(calFile), calFile);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }

    cdfStatus = CDFcreateAttr(cdf, "TEXT", GLOBAL_SCOPE, &attrNum);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    char *text = "CHAOS-7-based THEMIS ASI field-line trace (https://github.com/JohnathanBurchill/chaos).";
    cdfStatus = CDFputAttrgEntry(cdf, attrNum, entry, CDF_CHAR, strlen(text), text);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }

    // Variable attributes
    cdfStatus = CDFcreateAttr(cdf, "Name", VARIABLE_SCOPE, &attrNum);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    cdfStatus = CDFcreateAttr(cdf, "Description", VARIABLE_SCOPE, &attrNum);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    cdfStatus = CDFcreateAttr(cdf, "Unit", VARIABLE_SCOPE, &attrNum);
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }

    cdfStatus = addVariableAttributes(cdf, "SiteAbbreviation", "Camera site abbreviation", "-");
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    cdfStatus = addVariableAttributes(cdf, "SiteLatitude", "Camera site geodetic latitude", "degree");
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    cdfStatus = addVariableAttributes(cdf, "SiteLongitude", "Camera site longitude", "degree");
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    cdfStatus = addVariableAttributes(cdf, "SiteAltitude", "Camera site altitude", "m");
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    cdfStatus = addVariableAttributes(cdf, "EmissionLatitudes", "Geocentric latitudes of the THEMIS ASI pixel corners calculated for the assumed emission altitude.", "degree");
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    cdfStatus = addVariableAttributes(cdf, "EmissionLongitudes", "Geocentric longitudes of the THEMIS ASI pixel corners calculated for the assumed emission altitude.", "degree");
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    cdfStatus = addVariableAttributes(cdf, "EmissionRadii", "Geocentric radii of the THEMIS ASI pixel corners calculated for the assumed emission altitude.", "m");
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    cdfStatus = addVariableAttributes(cdf, "EmissionAltitudeGeodetic", "Assumed geodetic altitude of the auroral emissions.", "m");
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    cdfStatus = addVariableAttributes(cdf, "TargetLatitudes", "Geocentric latitudes of the THEMIS ASI pixel corners calculated by tracing along the CHAOS-7 magnetic lines of force to the target altitude.", "degree");
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    cdfStatus = addVariableAttributes(cdf, "TargetLongitudes", "Geocentric longitudes of the THEMIS ASI pixel corners calculated by tracing along the CHAOS-7 magnetic lines of force to the target altitude.", "degree");
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    cdfStatus = addVariableAttributes(cdf, "TargetRadii", "Geocentric radii of the THEMIS ASI pixel corners calculated by tracing along the CHAOS-7 magnetic lines of force to the target altitude.", "m");
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    cdfStatus = addVariableAttributes(cdf, "TargetAltitudeGeodetic", "Target altitude for field-line tracing.", "m");
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }
    cdfStatus = addVariableAttributes(cdf, "CalibrationEpoch", "Epoch of ASI calibration. CHAOS-7 is initialized to the beginning of the day.", "milliseconds from 0000-01-01:00:00:00 UT, no leap seconds");
    if (cdfStatus != CDF_OK)
    {
        CDFcloseCDF(cdf);
        return EXIT_FAILURE;
    }

    cdfStatus = CDFcloseCDF(cdf);
    if (cdfStatus != CDF_OK)
    {
        char msg[CDF_STATUSTEXT_LEN+1];
        CDFgetStatusText(cdfStatus, msg);
        fprintf(stderr, "Could not close exported CDF: %s\n", msg);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int addVariableAttributes(CDFid cdf, char *name, char *description, char *units)
{
    int cdfStatus = CDF_OK;
    long varNum = CDFgetVarNum(cdf, name);

    cdfStatus = CDFputAttrzEntry(cdf, CDFgetAttrNum(cdf, "Name"), varNum, CDF_CHAR, strlen(name), name);
    if (cdfStatus != CDF_OK)
        goto cleanup;

    cdfStatus = CDFputAttrzEntry(cdf, CDFgetAttrNum(cdf, "Description"), varNum, CDF_CHAR, strlen(description), description);
    if (cdfStatus != CDF_OK)
        goto cleanup;

    cdfStatus = CDFputAttrzEntry(cdf, CDFgetAttrNum(cdf, "Unit"), varNum, CDF_CHAR, strlen(units), units);
    if (cdfStatus != CDF_OK)
        goto cleanup;

cleanup:
    return cdfStatus;
}
