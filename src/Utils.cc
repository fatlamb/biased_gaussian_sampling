#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>

#include <iostream>

#include "Utils.h"

void print_gslvec(gsl_vector* vec)
{
	int i;
	for(i=0; i<vec->size; i++)
	{
		std::cout << gsl_vector_get(vec,i) << std::endl;
	}
}

void print_gslmat(gsl_matrix* mat)
{
	std::cout <<"#--------------------------------------------------#"<<std::endl;
	int i,j;
	for(i=0; i<mat->size1; i++)
	{
		for(j=0; j<mat->size2; j++)
		{
			std::cout << gsl_matrix_get(mat,i,j) << "  ";
		}

		std::cout << std::endl;
	}
	std::cout <<"#--------------------------------------------------#"<<std::endl;
}

void linear_vector_ramp(gsl_vector * vec, double r_final)
{
	int nsamples=vec->size;

	if (nsamples==1)
	{
		gsl_vector_set(vec,0,0);
	}
	else
	{
		int i;
		double rfrac;
		for(i=0; i<nsamples; i++)
		{
			rfrac=((double)i/((double)nsamples - 1));	
			gsl_vector_set(vec,i,r_final*rfrac);
		}
	}
}

