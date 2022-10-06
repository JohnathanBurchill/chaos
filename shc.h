/*

    CHAOS: shc.c

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

#ifndef _CHAOS_SHC_H
#define _CHAOS_SHC_H

#include <stdio.h>
#include <stdbool.h>

#define SHC_INFO_BUFFER_SIZE 1024

enum SHCError {
    SHC_OK = 0,
    SHC_DIRECTORY_READ,
    SHC_FILE_READ,
    SHC_FILE_CONTENTS,
    SHC_NUMBER_OF_COEFFICIENTS,
    SHC_MEMORY,
    SHC_INTERPOLATION,
    SHC_FRACTIONAL_YEAR

};

typedef struct SHCCoefficients
{
    bool initialized;
    const char coeffFilename[FILENAME_MAX];
    const char info[SHC_INFO_BUFFER_SIZE];
    int minimumN;
    int maximumN;
    int numberOfTimes;
    size_t numberOfTerms; // For GSL storage
    size_t gCoeffs;
    size_t hCoeffs;
    int bSplineOrder;
    int bSplineSteps;
    double *times;
    double *gTimeSeries;
    double *hTimeSeries;
    double *gNow;
    double *hNow;
    double *polynomials;
    double *derivatives;
    double *aoverrpowers;
} SHCCoefficients;

typedef struct ChaosCoefficients
{
    bool initialized;
    SHCCoefficients core;
    SHCCoefficients coreExtrapolation;
    SHCCoefficients crust;
} ChaosCoefficients;


int loadModelCoefficients(const char *coeffDir, ChaosCoefficients *coeffs);
int loadSHCCoefficients(SHCCoefficients *coeffs);

void freeChaosCoefficients(ChaosCoefficients *coeffs);
void freeSHCCoefficients(SHCCoefficients *coeffs);

int interpolateSHCCoefficients(ChaosCoefficients *coeffs, int year, int month, int day);

int yearFraction(long year, long month, long day, double* fractionalYear);


#endif // _CHAOS_SHC_H
