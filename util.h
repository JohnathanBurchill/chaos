#ifndef _CHAOS_UTIL_H
#define _CHAOS_UTIL_H

void elazToEnu(double elevation, double azimuth, double *xyz);
void geodeticToXyz(double geodeticLatitude, double longitude, double heightM, double *xyz);
void geodeticToGeocentric(double geodeticLatitude, double longitude, double heightM, double *geocentricLatitude, double *geocentricRadius, double *xyz);
void xyzToGeodetic(double *xyz, double *geodeticLatitude, double *longitude, double *heightM);
void lookDirectionToPosition(double *sitePositionXyz, double *lookDirectionEnu, double targetAltKm, double *geocentricPosition);

#endif // _CHAOS_UTIL_H
