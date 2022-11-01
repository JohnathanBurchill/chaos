#include "trace.h"

#include "shc.h"
#include "model.h"

#include <stdio.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>

int initializeTracer(char *coeffDir, int year, int month, int day, ChaosCoefficients *coeffs)
{
    int status = CHAOS_MODEL_OK;

	status = loadModelCoefficients(coeffDir, coeffs);
	if (status != SHC_OK || !coeffs->initialized)
        return status;

	status = interpolateSHCCoefficients(coeffs, year, month, day);
	if (status != SHC_OK)
    {
        freeChaosCoefficients(coeffs);
        return status;
    }
}

int trace(ChaosCoefficients *coeffs, int startingDirection, double latitude, double longitude, double alt1km, double minAltkm, double maxAltkm, double *latitude2, double *longitude2, double *altitude2, long *stepsTaken)
{

    if (latitude2 == NULL || longitude2 == NULL || altitude2 == NULL)
        return CHAOS_TRACE_POINTER;

    int status = CHAOS_TRACE_OK;

    // NEC system
    double bField[3] = {0.0};

    double degrees = M_PI / 180.0;

    double theta = (90.0 - latitude) * degrees;
    double phi = longitude * degrees;
    double earthRadiuskm = 6371.2;
    double r = earthRadiuskm + alt1km;
    // Test calculation to see if we can calculate field without error
    status = internalFieldNEC(r, theta, phi, coeffs, bField);
    if (status != CHAOS_MODEL_OK)
        return status;

    double rMin = earthRadiuskm + minAltkm;
    double rMax = earthRadiuskm + maxAltkm;

    size_t maxSteps = 10000;

    size_t steps = 0;
    // starting cartesian position
    // initial velocity is 0.0
    double y[3] = {0.0};
    y[0] = r * sin(theta) * cos(phi);
    y[1] = r * sin(theta) * sin(phi);
    y[2] = r * cos(theta);
    double yOld[3] = {0.0};
    for (int i = 0; i < 3; i++)
        yOld[i] = y[i];

    TracingState state = {0};
    state.coeffs = coeffs;
    state.startingDirection = (double) startingDirection; // +1 is parallel to B
    state.currentDirection = state.startingDirection;
    state.speed = 10.0; // km/s

    gsl_odeiv2_system sys = {force, NULL, 3, &state};
    double h = 0.5;
    double threshold=1e-2;
    // Something something "variable-coefficient linear multistep Adams method in Nordsieck form"
    // From GSL doc on Ordinary Differential Equations
    // Faster than RKF45
    // apex.f (A.D. Richmond, August 1994: 
    //   ITRACE subroutine: "This uses the 4-point Adams formula after initialization.
    //                       First 7 iterations advance point by 3 steps.""
    // Same Adams?

    // Need to check accuracy of the results at some point.
    gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_msadams, h, threshold, 0.0);
   
    double t = 0.0;
    double ti = 0.0;
    double dtMax = 0.5;
    double tOld = 0.0;
    double rOld = r;
    while (steps < maxSteps && r < rMax && r >= rMin)
    {
        status = gsl_odeiv2_driver_apply(driver, &t, t + dtMax, y);
        if (status != GSL_SUCCESS)
        {
            gsl_odeiv2_driver_free(driver);
            return CHAOS_TRACE_GSL_ERROR;
        }
        r = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
        // printf("h: %.1lf\t\n", r-earthRadiuskm);
        if ((r - rMax) > 0.0000001 || (r - rMin) < -0.0000001)
        {
            // printf("Resetting...\n");
            t = tOld;
            r = rOld;
            for (int i = 0; i < 3; i++)
                y[i] = yOld[i];
            h /= 2.0;
            dtMax /= 2.0;
            gsl_odeiv2_driver_reset(driver);
        }
        else
        {
            tOld = t;
            rOld = r;
            for (int i = 0; i < 3; i++)
                yOld[i] = y[i];
            steps++;
        }
    }

    // Store result
    theta = acos(y[2] / r);
    phi = atan2(y[1], y[0]);
    *latitude2 = 90.0 - theta / degrees;
    *longitude2 = phi / degrees;
    *altitude2 = r - earthRadiuskm;

    if (stepsTaken != NULL)
        *stepsTaken = steps;

    gsl_odeiv2_driver_free(driver);

    return CHAOS_TRACE_OK;

}

int force(double t, const double y[], double f[], void *data)
{
    (void)t;
    TracingState *s = (TracingState*)data;

    double r = 0.0;
    double theta = 0.0;
    double phi = 0.0;
    double startingDirection = s->startingDirection;
    double currentDirection = s->currentDirection;

    // BNEC
    double b[3] = {0.0};
    // BXYZ
    double bxyz[3] = {0.0};

    r = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
    theta = acos(y[2] / r);
    phi = atan2(y[1], y[0]);
    internalFieldNEC(r, theta, phi, s->coeffs, b);

    double n[3] = {0.0};
    double e[3] = {0.0};
    double c[3] = {0.0};

    c[0] = -sin(theta) * cos(phi);
    c[1] = -sin(theta) * sin(phi);
    c[2] = -cos(theta);

    e[0] = c[1] * 1.0;
    e[1] = -c[0] * 1.0;
    e[2] = 0.0;

    n[0] = e[1] * c[2] - e[2] * c[1];
    n[1] = -e[0] * c[2] + e[2] * c[0];
    n[2] = e[0] * c[1] - e[1] * c[0];

    double nMag = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    n[0] /= nMag;
    n[1] /= nMag;
    n[2] /= nMag;

    // Calculate bx, by, bz
    bxyz[0] = b[0] * n[0] + b[1] * e[0] + b[2] * c[0];       
    bxyz[1] = b[0] * n[1] + b[1] * e[1] + b[2] * c[1];       
    bxyz[2] = b[0] * n[2] + b[1] * e[2] + b[2] * c[2];       

    double bMag = sqrt(bxyz[0] * bxyz[0] + bxyz[1] * bxyz[1] + bxyz[2] * bxyz[2]);

    // Velocity is set to 1 km/s parallel to field (for currentDirection == 1)
    f[0] = s->speed * s->currentDirection * bxyz[0] / bMag;
    f[1] = s->speed * s->currentDirection * bxyz[1] / bMag;
    f[2] = s->speed * s->currentDirection * bxyz[2] / bMag;

    return GSL_SUCCESS;

}

int internalFieldNEC(double r, double theta, double phi, ChaosCoefficients *coeffs, double *bInt)
{
    double b[3] = {0.0, 0.0, 0.0};

    int status = CHAOS_MODEL_OK;
    status = calculateField(r, theta, phi, &coeffs->core, b, b+1, b+2);
    if (status != CHAOS_MODEL_OK)
        return status;

    *bInt = *b;
    *(bInt + 1) = *(b + 1);
    *(bInt + 2) = *(b + 2);

    status = calculateField(r, theta, phi, &coeffs->crust, b, b+1, b+2);
    if (status != CHAOS_MODEL_OK)
        return status;

    *bInt += *b;
    *(bInt + 1) += *(b + 1);
    *(bInt + 2) += *(b + 2);

    return CHAOS_MODEL_OK;
}