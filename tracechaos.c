#include "trace.h"

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
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

    if (argc - nOptions != 12)
    {
        printf("Incorrect number of arguments.\n");
        printf("usage: %s coeffDir tracingDirection year month day glat glon galt minAlt maxAlt altitudestep\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    char *coeffDir = argv[1];
    int startingDirection = atoi(argv[2]);
    int year = atoi(argv[3]);
    int month = atoi(argv[4]);
    int day = atoi(argv[5]);

    double latitude1 = atof(argv[6]);
    double longitude1 = atof(argv[7]);
    double altitudekm1 = atof(argv[8]);
    double minAlt = atof(argv[9]);
    double maxAlt = atof(argv[10]);
    double deltaAltkm = atof(argv[11]);

    double latitude2 = 0.0;
    double longitude2 = 0.0;

	ChaosCoefficients coeffs = {0};

    status = initializeTracer(coeffDir, year, month, day, &coeffs);


    // Does not work if we go too far along the field line. 
    // Good for high latitudes up to Swarm altitude
    long steps = 0;
    double finalAltitude = 0.0;
    for (double alt = altitudekm1; alt <= maxAlt; alt+=deltaAltkm)
    {
        // Inefficient. Could store steps along the way in the trace function
        status = trace(&coeffs, startingDirection, latitude1, longitude1, altitudekm1, minAlt, alt, &latitude2, &longitude2, &finalAltitude, &steps);

        printf("%lf %lf %lf %lf %lf %lf %ld\n", latitude1, longitude1, altitudekm1, latitude2, longitude2, finalAltitude, steps);

    }

    freeChaosCoefficients(&coeffs);

    return EXIT_SUCCESS;
}