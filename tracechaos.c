#include "trace.h"

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <gsl/gsl_errno.h>

sig_atomic_t keep_running;
char infoHeader[50];

int main (int argc, char *argv[])
{

    gsl_set_error_handler_off();
    sprintf(infoHeader, "%s", "tracechaos: ");
    keep_running = true;

    int status = CHAOS_TRACE_OK;

    int nOptions = 0;

    double minimumAltitudekm = 0.0;
    double accuracy = 0.001;

    for (int i = 0; i < argc; i++)
    {
        if (strncmp("--minimum-altitude-km=", argv[i], 22) == 0)
        {
            char *lastParsedChar = argv[i]+22;
            double value = strtod(argv[i] + 22, &lastParsedChar);
            if (lastParsedChar == argv[i] + 22)
            {
                fprintf(stderr, "%s: unable to parse %s\n", argv[0], argv[i]);
                exit(EXIT_FAILURE);
            }
            minimumAltitudekm = value;
            nOptions++;
        }
        if (strncmp("--accuracy=", argv[i], 11) == 0)
        {
            char *lastParsedChar = argv[i]+11;
            double value = strtod(argv[i] + 11, &lastParsedChar);
            if (lastParsedChar == argv[i] + 11)
            {
                fprintf(stderr, "%s: unable to parse %s\n", argv[0], argv[i]);
                exit(EXIT_FAILURE);
            }
            accuracy = value;
            nOptions++;
        }
    }


    if (argc - nOptions != 12)
    {
        printf("Incorrect number of arguments.\n");
        printf("usage: %s coeffDir tracingDirection year month day glat glon startAlt stopAlt1 stopAlt2 altitudeStep [--minimum-altitude-km=value]\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    char *coeffDir = argv[1];
    int startingDirection = atoi(argv[2]);
    int year = atoi(argv[3]);
    int month = atoi(argv[4]);
    int day = atoi(argv[5]);

    double latitude1 = atof(argv[6]);
    double longitude1 = atof(argv[7]);
    double startAlt = atof(argv[8]);
    double stopAlt1 = atof(argv[9]);
    double stopAlt2 = atof(argv[10]);
    double deltaAltkm = atof(argv[11]);

    double latitude2 = 0.0;
    double longitude2 = 0.0;

	ChaosCoefficients coeffs = {0};

    status = initializeTracer(coeffDir, year, month, day, &coeffs);


    // Does not work if we go too far along the field line. 
    // Good for high latitudes up to Swarm altitude
    long steps = 0;
    double finalAltitude = 0.0;
    for (double alt = stopAlt1; alt <= stopAlt2; alt+=deltaAltkm)
    {
        // Inefficient. Could store steps along the way in the trace function
        status = trace(&coeffs, startingDirection, accuracy, latitude1, longitude1, startAlt, minimumAltitudekm, alt, &latitude2, &longitude2, &finalAltitude, &steps);

        printf("%lf %lf %lf %lf %lf %lf %ld\n", latitude1, longitude1, startAlt, latitude2, longitude2, finalAltitude, steps);

    }

    freeChaosCoefficients(&coeffs);

    return EXIT_SUCCESS;
}