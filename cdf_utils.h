/*

    CHAOS: cdf_utils.h

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

int getOutputFilename(const char satellite, long year, long month, long day, const char *exportDir, double *beginTime, double *endTime, char *cdfFileName, char *magDataset);

CDFstatus exportCdf(const char *cdfFilename, const char *magFilename, const char *shcFile1, const char *shcFile2, const char satellite, const char *dataset, const char *exportVersion, double *times, double *latitudes, double *longitudes, double *radii, double *bCore, double *bCrust, double *dbMeas, size_t nVectors);

void exportMetaInfo(const char *outputFilename, const char *magFilename, const char *chaosCoreFilename, const char *chaosStaticFilename, long nVectors, time_t startTime, time_t stopTime);

enum EXPORT_FLAGS {
    EXPORT_OK = 0,
    EXPORT_MEM = 1
};


#endif // _LOAD_CDF_H

