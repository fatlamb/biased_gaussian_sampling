#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <iostream>

#include "RanGen.h"

void gen_normal_vector(gsl_vector* z)
{
	const gsl_rng_type * T;
	gsl_rng * r;

	int n=z->size;
	double sigma = 1.0;

	/* create a generator chosen by the 
	   environment variable GSL_RNG_TYPE */

	gsl_rng_env_setup();

	T = gsl_rng_default;
	r = gsl_rng_alloc (T);


//FIXME: use timestamp as seed??

//  printf ("generator type: %s\n", gsl_rng_name (r));
//  printf ("seed = %lu\n", gsl_rng_default_seed);
//  printf ("first value = %lu\n", gsl_rng_get (r));


	/* print n random variates chosen from 
	   the poisson distribution with mean 
	   parameter mu */

	int i;
	double el;
	for (i = 0; i < n; i++) 
	{
		el = gsl_ran_gaussian_ziggurat(r, sigma);
		gsl_vector_set(z,i,el);	
	}

	gsl_rng_free (r);
}
