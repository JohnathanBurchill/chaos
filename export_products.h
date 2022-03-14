#ifndef _EXPORT_PRODUCTS_H
#define _EXPORT_PRODUCTS_H

#include <stdint.h>
#include <time.h>

#include <cdf.h>

CDFstatus exportProducts(const char *slidemFilename, char satellite, double beginTime, double endTime, uint8_t **hmDataBuffers, long nHmRecs, double *vn, double *ve, double *vc, double *ionEffectiveMass, double *ionDensity, double *ionDriftRaw, double *ionDrift, double *ionEffectiveMassError, double *ionDensityError, double *ionDriftError, double *fpAreaOML, double *rProbeOML, double *electronTemperature, double *spacecraftPotential, double *ionEffectiveMassTTS, uint32_t *mieffFlags, uint32_t *viFlags, uint32_t *niFlags);

CDFstatus exportSlidemCdf(const char *cdfFilename, const char satellite, const char *exportVersion, uint8_t **hmDataBuffers, long nHmRecs, double *vn, double *ve, double *vc, double *ionEffectiveMass, double *ionDensity, double *ionDriftRaw, double *ionDrift, double *ionEffectiveMassError, double *ionDensityError, double *ionDriftError, double *fpAreaOML, double *rProbeOML, double *electronTemperature, double *spacecraftPotential, double *ionEffectiveMassTTS, uint32_t *mieffFlags, uint32_t *viFlags, uint32_t *niFlags);

void exportSlidemMetainfo(const char *slidemFilename, const char *fpFilename, const char *hmFilename, const char *magFilename, const char *modFilename, const char *modFilenamePrevious, long nVnecRecsPrev, time_t startTime, time_t stopTime);


enum EXPORT_FLAGS {
    EXPORT_OK = 0,
    EXPORT_MEM = 1
};

#endif // _EXPORT_PRODUCTS_H
