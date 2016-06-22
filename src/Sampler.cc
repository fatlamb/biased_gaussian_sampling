#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>

#include <iostream>

#include "Covgen.h"
#include "Utils.h"

int main()
{
	
	int nsamples = 10; //Total number, including origin sample.
	double r_final=10; //Final radius.  
	double k_cutoff=10; //Cutoff in spatial frequency.
	
	//Initialize sample radius vector.
	gsl_vector* sample_radii = gsl_vector_alloc(nsamples);	

	//Fill sample radius vector.
	linear_vector_ramp(sample_radii,r_final);
	

	//Initialize covariance matrix:

	gsl_matrix* cov_mat = gsl_matrix_alloc(nsamples,nsamples);
	
	print_gslmat(cov_mat);
	print_gslvec(sample_radii);

	calculate_cov(cov_mat,sample_radii,k_cutoff);
	print_gslmat(cov_mat);

	gsl_vector_free(sample_radii);
	gsl_matrix_free(cov_mat);
	return 0;
}
