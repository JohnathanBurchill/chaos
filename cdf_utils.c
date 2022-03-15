
#include "cdf_utils.h"
#include "chaos_settings.h"
#include "cdf_vars.h"

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <fts.h>
#include <cdf.h>
#include <signal.h>
#include <time.h>

extern volatile sig_atomic_t keep_running;

extern char infoHeader[50];

void loadCdf(const char *cdfFile, char *variables[], int nVariables, uint8_t **dataBuffers, size_t *numberOfRecords)
{
    // Open the CDF file with validation
    CDFsetValidate(VALIDATEFILEoff);
    CDFid cdfId;
    CDFstatus status;
    // Attributes
    long attrN;
    long entryN;
    char attrName[CDF_ATTR_NAME_LEN256+1];
    long attrScope, maxEntry;

    // Check CDF info
    long decoding, encoding, majority, maxrRec, numrVars, maxzRec, numzVars, numAttrs, format;

    long numBytesToAdd, numVarBytes, numValues;
    long varNum, dataType, numElems, numRecs, numDims, recVary;
    long dimSizes[CDF_MAX_DIMS], dimVarys[CDF_MAX_DIMS];
    CDFdata data;

    status = CDFopenCDF(cdfFile, &cdfId);
    if (status != CDF_OK) 
    {
        printErrorMessage(status);
        fprintf(stdout, "%s Could not open CDF file. Skipping this date.\n", infoHeader);
        return;
    }

    status = CDFgetFormat(cdfId, &format);
    status = CDFgetDecoding(cdfId, &decoding);
    status = CDFinquireCDF(cdfId, &numDims, dimSizes, &encoding, &majority, &maxrRec, &numrVars, &maxzRec, &numzVars, &numAttrs);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        fprintf(stdout, "\n%s Problem with CDF file. Skipping this date.\n", infoHeader);
        closeCdf(cdfId);
        return;
    }
    uint8_t nVars = numzVars;
    for (uint8_t i = 0; i<nVariables; i++)
    {
        status = CDFconfirmzVarExistence(cdfId, variables[i]);
        if (status != CDF_OK)
        {
            printErrorMessage(status);
            fprintf(stdout, "\n%s Error reading variable %s from CDF file. Skipping this date.\n", infoHeader, variables[i]);
            closeCdf(cdfId);
                return;
        }
    }
    

    for (uint8_t i = 0; i < nVariables && keep_running == 1; i++)
    {
        varNum = CDFgetVarNum(cdfId, variables[i]);
        status = CDFreadzVarAllByVarID(cdfId, varNum, &numRecs, &dataType, &numElems, &numDims, dimSizes, &recVary, dimVarys, &data);
        if (status != CDF_OK)
        {
            printErrorMessage(status);
            fprintf(stdout, "%s Error loading data for %s. Skipping this date.\n", infoHeader, variables[i]);
            closeCdf(cdfId);
            CDFdataFree(data);
            return;
        }
        // Calculate new size of memory to allocate
        status = CDFgetDataTypeSize(dataType, &numVarBytes);
        numValues = 1;
        for (uint8_t j = 0; j < numDims; j++)
        {
            numValues *= dimSizes[j];
        }
        numBytesToAdd = numValues * numRecs * numVarBytes;
        dataBuffers[i] = (uint8_t*) realloc(dataBuffers[i], (size_t) numBytesToAdd);
        memcpy(dataBuffers[i], data, numBytesToAdd);
        CDFdataFree(data);
    }
    // close CDF
    closeCdf(cdfId);

    // Update number of records found and memory allocated
    *numberOfRecords += numRecs;
    // *totalMemoryAllocated = fpMemorySize;

}

void printErrorMessage(CDFstatus status)
{
    char errorMessage[CDF_STATUSTEXT_LEN + 1];
    CDFgetStatusText(status, errorMessage);
    fprintf(stdout, "%s%s\n", infoHeader, errorMessage);
}

void closeCdf(CDFid id)
{
    CDFstatus status;
    status = CDFcloseCDF(id);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
    }
}

int getInputFilename(const char satelliteLetter, long year, long month, long day, const char *path, const char *dataset, char *filename)
{
    char *searchPath[2] = {NULL, NULL};
    searchPath[0] = (char *)path;

    FTS * fts = fts_open(searchPath, FTS_PHYSICAL | FTS_NOCHDIR, NULL);     
    if (fts == NULL)
    {
            printf("Could not open directory %s for reading.", path);
            return CDF_FIND_FILENAME;
    }
    FTSENT * f = fts_read(fts);

    bool gotFile = false;
    long fileYear;
    long fileMonth;
    long fileDay;
    long lastVersion = -1;
    long fileVersion;
    while(f != NULL)
    {
        // Most Swarm CDF file names have a length of 59 characters. The MDR_MAG_HR files have a length of 70 characters.
        // The MDR_MAG_HR files have the same filename structure up to character 55.
        if ((strlen(f->fts_name) == 59 || strlen(f->fts_name) == 70) && *(f->fts_name+11) == satelliteLetter && strncmp(f->fts_name+13, dataset, 5) == 0)
        {
            char fyear[5] = { 0 };
            char fmonth[3] = { 0 };
            char fday[3] = { 0 };
            char version[5] = { 0 };
            strncpy(fyear, f->fts_name + 19, 4);
            fileYear = atol(fyear);
            strncpy(fmonth, f->fts_name + 23, 2);
            fileMonth = atol(fmonth);
            strncpy(fday, f->fts_name + 25, 2);
            fileDay = atol(fday);
            strncpy(version, f->fts_name + 51, 4);
            fileVersion = atol(version);
            if (fileYear == year && fileMonth == month && fileDay == day && fileVersion > lastVersion)
            {
                lastVersion = fileVersion;
                sprintf(filename, "%s", f->fts_path);
                gotFile = true;
            }
        }
        f = fts_read(fts);
    }

    fts_close(fts);
    if (gotFile)
    {
        return CDF_ALL_GOOD;
    }
    else
    {
        return CDF_FIND_FILENAME;
    }

}

int getOutputFilename(const char satellite, long year, long month, long day, const char *exportDir, double *beginTime, double *endTime, char *cdfFileName)
{

    if (satellite != 'A' && satellite != 'B' && satellite != 'C')
    {
        return -1;
    }

    *beginTime = computeEPOCH(year, month, day, 0, 0, 0, 0);
    *endTime = computeEPOCH(year, month, day, 23, 59, 59, 999);

    sprintf(cdfFileName, "%s/SW_%s_MAG%c%s_2__%04d%02d%02dT000000_%04d%02d%02dT235959_%s", exportDir, CHAOS_PRODUCT_TYPE, satellite, CHAOS_PRODUCT_CODE, (int)year, (int)month, (int)day, (int)year, (int)month, (int)day, EXPORT_VERSION_STRING);

    return 0;

}


CDFstatus exportCdf(const char *cdfFilename, const char satellite, const char *exportVersion, double *times, double *bCore, double *bCrust, double *dbMeas, size_t nVectors)
{

    fprintf(stdout, "%sExporting CHAOS model data.\n",infoHeader);

    CDFid exportCdfId;
    CDFstatus status = CDF_OK;
    status = CDFcreateCDF((char *)cdfFilename, &exportCdfId);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        goto cleanup;
    }
    else
    {

        // export fpVariables
        createVarFrom1DVar(exportCdfId, "Timestamp", CDF_EPOCH, 0, nVectors-1, times);
        createVarFrom2DVar(exportCdfId, "B_core_nec", CDF_REAL8, 0, nVectors-1, bCore, 3);
        createVarFrom2DVar(exportCdfId, "B_crust_nec", CDF_REAL8, 0, nVectors-1, bCrust, 3);
        createVarFrom2DVar(exportCdfId, "dB_nec", CDF_REAL8, 0, nVectors-1, dbMeas, 3);

        // addAttributes(exportCdfId, SOFTWARE_VERSION_STRING, satellite, exportVersion, minTime, maxTime);

        fprintf(stdout, "%sExported %ld records to %s.cdf\n", infoHeader, nVectors, cdfFilename);
        fflush(stdout);
        status = CDF_OK;

    }

cleanup:
    closeCdf(exportCdfId);
    return status;

}


void exportMetaInfo(const char *outputFilename, const char *magFilename, const char *chaosCoreFilename, const char *chaosStaticFilename, long nVectors, time_t startTime, time_t stopTime)
{
    // Level 2 product ZIP file neads a HDR file, which is constructed from a metainfo file.
    char metaInfoFilename[FILENAME_MAX];
    sprintf(metaInfoFilename, "%s.metainfo", outputFilename);
    FILE *metaInfoFile = fopen(metaInfoFilename, "w");    
    if (metaInfoFile == NULL)
    {
        fprintf(stdout, "%sError opening metainfo file for writing.\n", infoHeader);
        return;
    }

    fprintf(metaInfoFile, "Type:%s\n", CHAOS_PRODUCT_TYPE);
    fprintf(metaInfoFile, "ProcessingCenter:UOC\n");
    fprintf(metaInfoFile, "Processor:UOC_CHAOS\n");
    fprintf(metaInfoFile, "ProcessorVersion:%s\n", SOFTWARE_VERSION);
    fprintf(metaInfoFile, "ProductError:0\n");

    fprintf(metaInfoFile, "Input:%s\n", magFilename + strlen(magFilename) - 70);
    fprintf(metaInfoFile, "Input:%s\n", chaosCoreFilename + strlen(chaosCoreFilename)-18);
    fprintf(metaInfoFile, "Input:%s\n", chaosStaticFilename + strlen(chaosStaticFilename)-20);

    char start[20] = {0};
    struct tm * tstart = gmtime(&startTime);
    sprintf(start, "%4d-%02d-%02dT%02d:%02d:%02d", tstart->tm_year+1900, tstart->tm_mon+1, tstart->tm_mday, tstart->tm_hour, tstart->tm_min, tstart->tm_sec);
    fprintf(metaInfoFile, "ProcessStart:%s\n", start);
    
    char stop[20] = {0};
    struct tm * tstop = gmtime(&stopTime);
    sprintf(stop, "%4d-%02d-%02dT%02d:%02d:%02d", tstop->tm_year+1900, tstop->tm_mon+1, tstop->tm_mday, tstop->tm_hour, tstop->tm_min, tstop->tm_sec);
    fprintf(metaInfoFile, "ProcessStop:%s\n", stop);

    size_t sLen = strlen(outputFilename);
    fprintf(metaInfoFile, "Output:%s.cdf\n", outputFilename + strlen(outputFilename) - 55);

    fclose(metaInfoFile);
    fprintf(stdout, "%sMetainfo file: %s\n", infoHeader, metaInfoFilename);

    return;


}