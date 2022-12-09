/*

    CHAOS: util.c

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

#include "util.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void elazToEnu(double elevation, double azimuth, double *enu)
{
    if (enu == NULL)
        return;

    double degree = M_PI / 180.0;
    double el = elevation * degree;
    double phi = fmod(90.0 - azimuth, 360.0) * degree;
    double ce = cos(el);

    enu[0] = cos(phi) * ce;
    enu[1] = sin(phi) * ce;
    enu[2] = sin(el);

    return;
}

void geodeticToGeocentric(double geodeticLatitude, double longitude, double heightM, double *geocentricLatitude, double *geocentricRadius, double *xyz)
{
    double a = 6378137.0;
    double b = 6356752.31425;
    double a2 = a * a;
    double b2 = b * b;

    double e2 = (a2 - b2) / a2;

    double sl = sin(geodeticLatitude * M_PI / 180.0);
    double cl = cos(geodeticLatitude * M_PI / 180.0);
    double cp = cos(longitude * M_PI / 180.0);
    double sp = sin(longitude * M_PI / 180.0);

    double sl2 = sl * sl;

    double n = a / sqrt(1.0 - e2 * sl2);

    double f = (a - b) / a;

    double nf = (n * (1 - f) * (1 - f) + heightM);
    double nh = n + heightM;

    double x = nh * cl * cp;
    double y = nh * cl * sp;
    double z = nf * sl;

    if (xyz != NULL)
    {
        xyz[0] = x;
        xyz[1] = y;
        xyz[2] = z;
    }

    double r = sqrt(x * x + y * y + z * z);

    double geocentricLat = atan(nf / nh * tan(geodeticLatitude * M_PI / 180.0)) * 180.0 / M_PI;

    if (geocentricLatitude != NULL)
        *geocentricLatitude = geocentricLat;

    if (geocentricRadius != NULL)
        *geocentricRadius = r;

    return;
}


void xyzToGeodetic(double *xyz, double *geodeticLatitude, double *longitude, double *heightM)
{
    if (xyz == NULL)
        return;

    double degree = M_PI / 180.0;

    if (longitude != NULL)
        *longitude = atan2(xyz[1], xyz[0]) / degree;

    // Newton Raphson iteration, https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_ECEF_to_geodetic_coordinates


    double a = 6378137.0;
    double b = 6356752.31425;
    double a2 = a * a;
    double b2 = b * b;
    double e2 = (a2 - b2) / a2;
    double kappa = 1.0 / (1 - e2);
    double kappa0 = kappa;
    double kappaOld = 0.0;

    double maxError = 0.0000000001;

    int nIterations = 0;
    double p = sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1]);
    double z = xyz[2];
    double c0 = 0.0;
    double c1 = 0.0;
    double p2 = p * p;
    double z2 = z * z;
    double kappa2 = kappa * kappa;
    while (fabs(kappa - kappaOld) > maxError)
    {
        nIterations++;
        kappaOld = kappa;
        p2 = p * p;
        z2 = z * z;
        kappa2 = kappa * kappa;
        c0 = p2 + (1-e2) * z2 * kappa2;
        c1 = sqrt(c0 * c0 * c0) / a / e2;
        kappa = 1.0 + (p2 + (1.0 - e2)*z2 * kappa2*kappa) / (c1 - p2);
        // printf("kappa0: %.9lf, kappa: %.9lf, kappa - kappaOld: %.10lf, c1: %lf\n", kappa0, kappa, kappa-kappaOld, c1);
    }
    // printf("nIter: %d\n", nIterations);

    if (geodeticLatitude != NULL)
        *geodeticLatitude = atan(kappa * z / p) / degree;

    if (heightM != NULL)
        *heightM = (1.0 / kappa - 1.0 / kappa0) * sqrt(p2 + z2*kappa2) / e2;

    return;
}

void lookDirectionToPosition(double *sitePositionXyz, double *lookDirectionEnu, double targetAltKm, double *geocentricPosition)
{
    if (sitePositionXyz == NULL || lookDirectionEnu == NULL)
        return;

    if (lookDirectionEnu[2] <= 0.0)
    {
        // Return NaNs for look directions below the horizon
        if (geocentricPosition != NULL)
        {
            geocentricPosition[0] = nan("");
            geocentricPosition[1] = nan("");
            geocentricPosition[2] = nan("");
        }
        return;
    }

    double rEarth = 6371200.0;
    double targetAltM = 1000.0 * targetAltKm;

    double length = sqrt(sitePositionXyz[0] * sitePositionXyz[0] + sitePositionXyz[1] * sitePositionXyz[1] + sitePositionXyz[2] * sitePositionXyz[2]);

    double uhat[3] = {0};
    for (int k = 0; k < 3; k++)
        uhat[k] = sitePositionXyz[k] / length;

    double ehat[3] = {0};
    ehat[0] = -uhat[1];
    ehat[1] = uhat[0];
    ehat[2] = 0.0;
    length = sqrt(ehat[0] * ehat[0] + ehat[1] * ehat[1] + ehat[2] * ehat[2]);
    for (int k = 0; k < 3; k++)
        ehat[k] /= length;

    double nhat[3] = {0};
    nhat[0] = uhat[1] * ehat[2] - uhat[2] * ehat[1];
    nhat[1] = -uhat[0] * ehat[2] + uhat[2] * ehat[0];
    nhat[2] = uhat[0] * ehat[1] - uhat[1] * ehat[0];
    length = sqrt(nhat[0] * nhat[0] + nhat[1] * nhat[1] + nhat[2] * nhat[2]);
    for (int k = 0; k < 3; k++)
        nhat[k] /= length;

    double el = atan(lookDirectionEnu[2] / sqrt(lookDirectionEnu[0]*lookDirectionEnu[0] + lookDirectionEnu[1] * lookDirectionEnu[1]));

    // Project a ray along the look direction by about the right distance
    // Calculate xyz position of ray 
    // Calculate geodetic altitude
    // Iterate by updating distance along the look direction 
    // until geodetic altitude == targetAltitude

    // Approximate distance is given by d^2 + 2 R d sin(el) == h^2 + 2 h R
    // with R - Earth radius (assuming a spherical Earth)
    //      h - altitude
    //     el - elevation angle
    float b = 2.0 * rEarth * sin(el);
    float c = - (targetAltM * targetAltM + 2.0 * targetAltM * rEarth);
    double projectionDistanceM = (-b + sqrt(b * b - 4.0 * c)) / 2.0;
    double posXyz[3] = {0};
    double geodeticAltM = targetAltM - 1000.0;

    // threshold
    double maxAltErrorM = 1.0;
    int nIterations = 0;
    double distanceStepM = 5000.0;
    double longitude = 0.0;
    double geodeticLatitude = 0.0;
    int lastSign = -1;
    int sign = -1;
    while (abs(geodeticAltM - targetAltM) > maxAltErrorM)
    {
        sign = geodeticAltM >= targetAltM ? 1 : -1;
        // short
        if (sign != lastSign)
            distanceStepM *= 0.25;
        lastSign = sign;
        if (geodeticAltM > targetAltM)
            projectionDistanceM -= distanceStepM;
        else
            projectionDistanceM += distanceStepM;

        for (int k = 0; k < 3; k++)
            posXyz[k] = sitePositionXyz[k] + projectionDistanceM * (lookDirectionEnu[0] * ehat[k] + lookDirectionEnu[1] * nhat[k] + lookDirectionEnu[2] * uhat[k]);

        xyzToGeodetic(posXyz, &geodeticLatitude, &longitude, &geodeticAltM);
        // printf("Approx dist: %.2lf, Geodetic alt: %.2lf, targetAlt: %.2lf, Geodetic alt - targetAlt: %.2lf\n", projectionDistanceM/1000.0, geodeticAltM/1000.0, targetAltM/1000.0, (geodeticAltM - targetAltM)/1000.0);
        nIterations++;

    }
    // printf("nIter: %d\n", nIterations);

    if (geocentricPosition != NULL)
    {
        // Lat, lon, radius
        geocentricPosition[0] = atan(posXyz[2] / sqrt(posXyz[0] * posXyz[0] + posXyz[1] * posXyz[1])) / M_PI * 180.0;
        geocentricPosition[1] = longitude;
        geocentricPosition[2] = sqrt(posXyz[0] * posXyz[0] + posXyz[1] * posXyz[1] + posXyz[2] * posXyz[2]);
    }

    return;
}
