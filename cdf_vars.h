#ifndef CDF_VARS_H 
#define CDF_VARS_H

#include <stdint.h>

#include <cdf.h>

CDFstatus createVarFrom1DVar(CDFid id, char *name, long dataType, long startIndex, long stopIndex, void *buffer);
CDFstatus createVarFrom2DVar(CDFid id, char *name, long dataType, long startIndex, long stopIndex, void *buffer1D, uint8_t dimSize);

#endif // CDF_VARS_H
