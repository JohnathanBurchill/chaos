#ifndef _LOAD_CDF_H
#define _LOAD_CDF_H

#include <stdint.h>
#include <stdlib.h>
#include <cdf.h>

enum CDF_UTILS {
    CDF_ALL_GOOD = 0,
    CDF_FIND_FILENAME = -1
};

void loadCdf(const char *cdfFile, char *variables[], int nVariables, uint8_t **dataBuffers, size_t *numberOfRecords);

void printErrorMessage(CDFstatus status);

void closeCdf(CDFid id);

int getInputFilename(const char satelliteLetter, long year, long month, long day, const char *path, const char *dataset, char *filename);

int getOutputFilename(const char satellite, long year, long month, long day, const char *exportDir, double *beginTime, double *endTime, char *cdfFileName);

CDFstatus exportCdf(const char *cdfFilename, const char satellite, const char *exportVersion, double *times, double *bCore, double *bCrust, double *bMeas, size_t nVectors);

void exportMetaInfo(const char *outputFilename, const char *magFilename, const char *chaosCoreFilename, const char *chaosStaticFilename, long nVectors, time_t startTime, time_t stopTime);

enum EXPORT_FLAGS {
    EXPORT_OK = 0,
    EXPORT_MEM = 1
};


#endif // _LOAD_CDF_H

