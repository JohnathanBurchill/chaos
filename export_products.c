/*

    CHAOS: export_products.c

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

#include "export_products.h"

#include "slidem.h"
#include "slidem_settings.h"
#include "cdf_vars.h"
#include "cdf_attrs.h"
#include "utilities.h"

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

// #include <zip.h>
#include <cdf.h>


extern char infoHeader[50];

CDFstatus exportProducts(const char *slidemFilename, char satellite, double beginTime, double endTime, uint8_t **hmDataBuffers, long nHmRecs, double *vn, double *ve, double *vc, double *ionEffectiveMass, double *ionDensity, double *ionDriftRaw, double *ionDrift, double *ionEffectiveMassError, double *ionDensityError, double *ionDriftError, double *fpAreaOML, double *rProbeOML, double *electronTemperature, double *spacecraftPotential, double *ionEffectiveMassTTS, uint32_t *mieffFlags, uint32_t *viFlags, uint32_t *niFlags)
{
    long hmTimeIndex = 0;
    beginTime = HMTIME();
    hmTimeIndex = nHmRecs - 1;
    endTime = HMTIME();

    CDFstatus status = CDF_OK;

    status = exportSlidemCdf(slidemFilename, satellite, EXPORT_VERSION_STRING, hmDataBuffers, nHmRecs, vn, ve, vc, ionEffectiveMass, ionDensity, ionDriftRaw, ionDrift, ionEffectiveMassError, ionDensityError, ionDriftError, fpAreaOML, rProbeOML, electronTemperature, spacecraftPotential, ionEffectiveMassTTS, mieffFlags, viFlags, niFlags);
    if (status != CDF_OK)
    {
        return status;
    }

    double minutesExported = (endTime - beginTime)/1000./60.;

    // report
    fprintf(stdout, "%sExported ~%.0f orbits (%ld 2 Hz records) of SLIDEM IDM data. %.1f%% coverage.\n", infoHeader, minutesExported/94., nHmRecs, minutesExported/1440.0*100.0); 

    return status;
}

CDFstatus exportSlidemCdf(const char *cdfFilename, const char satellite, const char *exportVersion, uint8_t **hmDataBuffers, long nHmRecs, double *vn, double *ve, double *vc, double *ionEffectiveMass, double *ionDensity, double *ionDriftRaw, double *ionDrift, double *ionEffectiveMassError, double *ionDensityError, double *ionDriftError, double *fpAreaOML, double *rProbeOML, double *electronTemperature, double *spacecraftPotential, double *ionEffectiveMassTTS, uint32_t *mieffFlags, uint32_t *viFlags, uint32_t *niFlags)
{

    fprintf(stdout, "%sExporting SLIDEM IDM data.\n",infoHeader);

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
        // velocity is a 1D variable (scalars are 0D in CDF parlance), per request of DTU
        double * vnec = malloc((size_t) (nHmRecs * 3 * sizeof(double)));
        if (vnec == NULL)
        {
            fprintf(stdout, "%s could not allocate memory to store VNEC.\n", infoHeader);
            status = EXPORT_MEM;
            goto cleanup;
        }
        for (long hmTimeIndex = 0; hmTimeIndex < nHmRecs; hmTimeIndex++)
        {
            vnec[3*hmTimeIndex] = vn[hmTimeIndex];
            vnec[3*hmTimeIndex + 1] = ve[hmTimeIndex];
            vnec[3*hmTimeIndex + 2] = vc[hmTimeIndex];
        }

        // export fpVariables
        createVarFrom1DVar(exportCdfId, "Timestamp", CDF_EPOCH, 0, nHmRecs-1, hmDataBuffers[0]);
        createVarFrom1DVar(exportCdfId, "Latitude", CDF_REAL8, 0, nHmRecs-1, hmDataBuffers[1]);
        createVarFrom1DVar(exportCdfId, "Longitude", CDF_REAL8, 0, nHmRecs-1, hmDataBuffers[2]);
        createVarFrom1DVar(exportCdfId, "Radius", CDF_REAL8, 0, nHmRecs-1, hmDataBuffers[3]);
        createVarFrom1DVar(exportCdfId, "Height", CDF_REAL8, 0, nHmRecs-1, hmDataBuffers[4]);
        createVarFrom1DVar(exportCdfId, "QDLatitude", CDF_REAL8, 0, nHmRecs-1, hmDataBuffers[5]);
        createVarFrom1DVar(exportCdfId, "MLT", CDF_REAL8, 0, nHmRecs-1, hmDataBuffers[7]);
        createVarFrom2DVar(exportCdfId, "V_sat_nec", CDF_REAL8, 0, nHmRecs-1, vnec, 3);
        createVarFrom1DVar(exportCdfId, "M_i_eff", CDF_REAL8, 0, nHmRecs-1, ionEffectiveMass);
        createVarFrom1DVar(exportCdfId, "M_i_eff_err", CDF_REAL8, 0, nHmRecs-1, ionEffectiveMassError);
        createVarFrom1DVar(exportCdfId, "M_i_eff_Flags", CDF_UINT4, 0, nHmRecs-1, mieffFlags);
        createVarFrom1DVar(exportCdfId, "M_i_eff_tbt_model", CDF_REAL8, 0, nHmRecs-1, ionEffectiveMassTTS);
        createVarFrom1DVar(exportCdfId, "V_i", CDF_REAL8, 0, nHmRecs-1, ionDrift);
        createVarFrom1DVar(exportCdfId, "V_i_err", CDF_REAL8, 0, nHmRecs-1, ionDriftError);
        createVarFrom1DVar(exportCdfId, "V_i_Flags", CDF_UINT4, 0, nHmRecs-1, viFlags);
        createVarFrom1DVar(exportCdfId, "V_i_raw", CDF_REAL8, 0, nHmRecs-1, ionDriftRaw);
        createVarFrom1DVar(exportCdfId, "N_i", CDF_REAL8, 0, nHmRecs-1, ionDensity);
        createVarFrom1DVar(exportCdfId, "N_i_err", CDF_REAL8, 0, nHmRecs-1, ionDensityError);
        createVarFrom1DVar(exportCdfId, "N_i_Flags", CDF_UINT4, 0, nHmRecs-1, niFlags);
        createVarFrom1DVar(exportCdfId, "A_fp", CDF_REAL8, 0, nHmRecs-1, fpAreaOML);
        createVarFrom1DVar(exportCdfId, "R_p", CDF_REAL8, 0, nHmRecs-1, rProbeOML);
        createVarFrom1DVar(exportCdfId, "T_e", CDF_REAL8, 0, nHmRecs-1, electronTemperature);
        createVarFrom1DVar(exportCdfId, "Phi_sc", CDF_REAL8, 0, nHmRecs-1, spacecraftPotential);

        // add attributes
        long hmTimeIndex = 0;
        double minTime = HMTIME();
        hmTimeIndex = nHmRecs - 1;
        double maxTime = HMTIME();

        addAttributes(exportCdfId, SOFTWARE_VERSION_STRING, satellite, exportVersion, minTime, maxTime);

        fprintf(stdout, "%sExported %ld records to %s.cdf\n", infoHeader, nHmRecs, cdfFilename);
        fflush(stdout);
        status = CDF_OK;

    }

cleanup:
    closeCdf(exportCdfId);
    return status;

}

void exportSlidemMetainfo(const char *slidemFilename, const char *fpFilename, const char *hmFilename, const char *magFilename, const char *modFilename, const char *modFilenamePrevious, long nVnecRecsPrev, time_t startTime, time_t stopTime)
{
    // Level 2 product ZIP file neads a HDR file, which is constructed from a metainfo file.
    char metaInfoFilename[FILENAME_MAX];
    sprintf(metaInfoFilename, "%s.metainfo", slidemFilename);
    FILE *metaInfoFile = fopen(metaInfoFilename, "w");    
    if (metaInfoFile == NULL)
    {
        fprintf(stdout, "%sError opening metainfo file for writing.\n", infoHeader);
        return;
    }

    fprintf(metaInfoFile, "Type:%s\n", SLIDEM_PRODUCT_TYPE);
    fprintf(metaInfoFile, "ProcessingCenter:UOC\n");
    fprintf(metaInfoFile, "Processor:UOC_SLIDEM\n");
    fprintf(metaInfoFile, "ProcessorVersion:%s\n", SOFTWARE_VERSION);
    fprintf(metaInfoFile, "ProductError:0\n");

    size_t fpLen = strlen(fpFilename);
    fprintf(metaInfoFile, "Input:%s\n", fpFilename + fpLen - 59);
    
    size_t hmLen = strlen(hmFilename);
    fprintf(metaInfoFile, "Input:%s\n", hmFilename + hmLen - 59);

    size_t vnecLen;
    if (nVnecRecsPrev > 0)
    {
        // used previous day MOD file to get coverage at beginning of processed day because
        // the MOD files have an offset of some seconds
        vnecLen = strlen(modFilenamePrevious);
        fprintf(metaInfoFile, "Input:%s\n", modFilenamePrevious + vnecLen - 59);
    }
    vnecLen = strlen(modFilename);
    fprintf(metaInfoFile, "Input:%s\n", modFilename + vnecLen - 59);

    // MAG CDF file
    vnecLen = strlen(magFilename);
    fprintf(metaInfoFile, "Input:%s\n", magFilename + vnecLen - 70);

    // F10.7
    fprintf(metaInfoFile, "Input:%s\n", "apf107.dat");

    // modified OML file
    fprintf(metaInfoFile, "Input:%s\n", ".slidem_modified_oml_configrc");

    char start[20] = {0};
    struct tm * tstart = gmtime(&startTime);
    sprintf(start, "%4d-%02d-%02dT%02d:%02d:%02d", tstart->tm_year+1900, tstart->tm_mon+1, tstart->tm_mday, tstart->tm_hour, tstart->tm_min, tstart->tm_sec);
    fprintf(metaInfoFile, "ProcessStart:%s\n", start);
    
    char stop[20] = {0};
    struct tm * tstop = gmtime(&stopTime);
    sprintf(stop, "%4d-%02d-%02dT%02d:%02d:%02d", tstop->tm_year+1900, tstop->tm_mon+1, tstop->tm_mday, tstop->tm_hour, tstop->tm_min, tstop->tm_sec);
    fprintf(metaInfoFile, "ProcessStop:%s\n", stop);

    size_t sLen = strlen(slidemFilename);
    fprintf(metaInfoFile, "Output:%s.cdf\n", slidemFilename + sLen - 55);

    fclose(metaInfoFile);
    fprintf(stdout, "%sMetainfo file: %s\n", infoHeader, metaInfoFilename);

    return;
}
