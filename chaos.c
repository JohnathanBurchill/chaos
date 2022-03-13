#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_math.h>


int main ()
{

	int status = 0;

	printf("Hi\n");

	double degrees = M_PI / 180.0;
	double a = 6371.2;
	double alt = 500.0;
	double r = a + alt;
	double theta = 67.0 * degrees;
	double phi = 31. * degrees;

	double *polynomials = NULL;
	double *aoverrpowers = NULL;
	double *times = NULL;
	double *gnm = NULL;
	double *hnm = NULL;

	char *coreFile = "shc/CHAOS-7.9_core.shc";
	char line[255];
	char character = '\0';

	FILE *f = fopen(coreFile, "r");	
	if (f == NULL)
	{
		printf("Could not open core model SHC file.\n");
		goto cleanup;
	}

	int minN, maxN, nTimes, bSplineOrder, bSplineSteps; 
	while (fgets(line, 255, f) != NULL && line[0] == '#');
	sscanf(line, "%d %d %d %d %d", &minN, &maxN, &nTimes, &bSplineOrder, &bSplineSteps);

	printf("MaxN: %d, nTimes: %d\n", maxN, nTimes);

	size_t n = gsl_sf_legendre_array_n(maxN);
	polynomials = malloc(n * sizeof(double));
	aoverrpowers = malloc(maxN * sizeof(double));
	times = (double*)malloc(nTimes * sizeof(double));
	// More memory than needed. Could be revised
	gnm = (double*)malloc(maxN * (maxN + 2) * nTimes * sizeof(double));
	hnm = (double*)malloc(maxN * (maxN + 2) * nTimes * sizeof(double));
	if (polynomials == NULL || aoverrpowers == NULL || times == NULL || gnm == NULL || hnm == NULL)
	{
		printf("Cannot remember anything. Check my memory.\n");
		goto cleanup;
	}

	int il, im;
	for (int i = 0; i < nTimes; i++)
	{
		fscanf(f, "%lf", times + i);
	}
	size_t gRead = 0;
	size_t hRead = 0;
	for (int i = 0; i < maxN * (maxN+2); i++)
	{
		fscanf(f, "%d %d", &il, &im);
		if (im >= 0)
		{
			for (int i = 0; i < nTimes; i++)
			{
				fscanf(f, "%lf", gnm + gRead);
				gRead++;
			}
		}
		else 
		{
			for (int i = 0; i < nTimes; i++)
			{
				fscanf(f, "%lf", hnm + hRead);
				hRead++;
			}
		}
	}
	// for (int i = 0; i < nTimes; i++)
	// {
	// 	printf("t = %lf\n", times[i]);
	// }
	// gRead = 0;
	// hRead = 0;
	// for (int l = 1; l <= maxN; l++)
	// {
	// 	for (int m = 0; m <= l; m++)
	// 	{
	// 		printf("g_%d_%d = %lf\n", l, m, gnm[gRead*nTimes + 0]);
	// 		gRead++;
	// 		if (m > 0)
	// 		{
	// 			printf("h_%d_%d = %lf\n", l, m, hnm[hRead*nTimes + 0]);
	// 			hRead++;
	// 		}
	// 	}
	// }
	fclose(f);
	f = NULL;

	double x = 0;
	double magneticPotential = 0.0;
	double magneticPotentialN = 0.0;

	// 2 Hz simulated
	double aoverr = 1.0;
	for (double t = 0.; t < 86400.; t+=0.5)
	{
		status = gsl_sf_legendre_array(GSL_SF_LEGENDRE_SCHMIDT, maxN, cos(theta), polynomials);
		if (status)
		{
			printf("status: %s\n", gsl_strerror(status));
			goto cleanup;
		}
		aoverr = a / r;
		aoverrpowers[0] = aoverr * aoverr; // (a/r)^n+1, n starting at 1
		for (int i = 2; i <= maxN; i++)
		{
			aoverrpowers[i-1] = aoverrpowers[i-2] * aoverr;
		}

		magneticPotential = 0.0;
		gRead = 0;
		hRead = 0;
		for (int l = 1; l <= maxN; l++)
		{
			magneticPotentialN = gnm[gRead*nTimes] * polynomials[gsl_sf_legendre_array_index(l, 0)];
			gRead++;
			for (int m = 1; m <= l; m++)
			{
				magneticPotentialN += (gnm[gRead*nTimes + 0] * cos(phi) + hnm[hRead*nTimes + 0] * sin(phi) ) * polynomials[gsl_sf_legendre_array_index(l, m)];
				gRead++;
				hRead++;
			}
			magneticPotentialN *= aoverrpowers[l-1] * a;
		}
		magneticPotential += magneticPotentialN;
	}
	// for (int i = 1; i <= maxN; i++)
	// 	printf("aoverr[%d]: %lf\n", i, aoverrpowers[i-1]);
	printf("magneticPotential(r=%lf, theta=%lf, phi=%lf) = %lf\n", r, theta, phi, magneticPotential);

cleanup:
	if (polynomials != NULL) free(polynomials);
	if (aoverrpowers != NULL) free(aoverrpowers);
	if (times != NULL) free(times);
	if (gnm != NULL) free(gnm);
	if (hnm != NULL) free(hnm);
	if (f != NULL) fclose(f);

	return 0;
}

