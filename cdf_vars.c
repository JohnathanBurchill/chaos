
#include "cdf_vars.h"
#include "cdf_utils.h"
#include "chaos_settings.h"

#include <stdlib.h>

#include <ctype.h>


CDFstatus createVarFrom1DVar(CDFid id, char *name, long dataType, long startIndex, long stopIndex, void *buffer)
{
    CDFstatus status;
    long exportDimSizes[1] = {0};
    long recVary = {VARY};
    long dimNoVary = {NOVARY};
    long varNumber;
    long cType;
    long cParams[CDF_MAX_PARMS];

    status = CDFcreatezVar(id, name, dataType, 1, 0L, exportDimSizes, recVary, dimNoVary, &varNumber);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        return status;
    }
    status = CDFsetzVarSparseRecords(id, varNumber, NO_SPARSERECORDS);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        return status;
    }
    cType = GZIP_COMPRESSION;
    cParams[0] = CDF_GZIP_COMPRESSION_LEVEL; // GZIP compression level 6 as suggested compromised between speed and size by CDF C reference 
    status = CDFsetzVarCompression(id, varNumber, cType, cParams);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        return status;
    }
    // Blocking factor 43200 as requested by DTU
    status = CDFsetzVarBlockingFactor(id, varNumber, CDF_BLOCKING_FACTOR);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        return status;
    }

    long dataTypeSize;
    status = CDFgetDataTypeSize(dataType, &dataTypeSize);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        return status;
    }
    status = CDFputVarRangeRecordsByVarName(id, name, 0, stopIndex-startIndex, (void*)((uint8_t*)buffer + (dataTypeSize*startIndex)));
    if (status != CDF_OK)
    {
        printErrorMessage(status);
    }
    return status;
}

CDFstatus createVarFrom2DVar(CDFid id, char *name, long dataType, long startIndex, long stopIndex, void *buffer1D, uint8_t dimSize)
{
    CDFstatus status = CDF_OK;
    long long nRecs = stopIndex - startIndex + 1;
    long dimSizes[1] = {0};
    long recVary = {VARY};
    long dimVary[1] = {VARY};
    long varNumber;
    long cType;
    long cParams[CDF_MAX_PARMS];

    dimSizes[0] = dimSize;

    status = CDFcreatezVar(id, name, dataType, 1, 1L, dimSizes, recVary, dimVary, &varNumber);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        return status;
    }
    cType = GZIP_COMPRESSION;
    cParams[0] = CDF_GZIP_COMPRESSION_LEVEL; // GZIP compression level 6 as suggested compromised between speed and size by CDF C reference 
    status = CDFsetzVarCompression(id, varNumber, cType, cParams);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        return status;
    }
    // Blocking factor 43200 as requested by DTU
    status = CDFsetzVarBlockingFactor(id, varNumber, CDF_BLOCKING_FACTOR);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        return status;
    }
    status = CDFsetzVarSparseRecords(id, varNumber, NO_SPARSERECORDS);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        return status;
    }
    long dataSize;
    status = CDFgetDataTypeSize(dataType, &dataSize);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
        return status;
    }
    status = CDFputVarRangeRecordsByVarName(id, name, 0, stopIndex-startIndex, (void *)buffer1D);
    if (status != CDF_OK)
    {
        printErrorMessage(status);
    }

    return status;
}

