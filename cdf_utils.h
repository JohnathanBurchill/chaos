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

#endif // _LOAD_CDF_H

