#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>


int main ()
{
	printf("Hi\n");

	size_t maxn = 20;

	double x = 0.0;

	size_t n = gsl_sf_legendre_array_n(maxn);

	double *polynomials = malloc(n * sizeof(double));

	int status = gsl_sf_legendre_array(GSL_SF_LEGENDRE_SCHMIDT, maxn, x, polynomials);

	printf("status: %s\n", gsl_strerror(status));

	for (size_t l = 1; l <= maxn; l++)
	{
		for (size_t m = 0; m <=l; m++)
		{
			printf("P_%ld_%ld(%f) = %lf\n", l, m, x, polynomials[gsl_sf_legendre_array_index(l, m)]);
		}
	}
	free(polynomials);

	return 0;
}

