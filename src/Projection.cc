#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <math.h>

#include <iostream>

#include "Projection.h"
#include "Utils.h"


void gslvector_outer_product(gsl_vector* vec1, gsl_vector* vec2, gsl_matrix* prod)
{
	int nrows=vec1->size;
	int ncols=vec2->size;

	double val1,val2;
	int i,j;
	for (i=0; i<nrows; i++)
	{
		val1=gsl_vector_get(vec1,i);	

		for (j=0; j<ncols; j++)
		{
			val2=gsl_vector_get(vec2,j);	
			gsl_matrix_set(prod,i,j,val1*val2);
		}
	}
}			

void construct_projector(gsl_matrix* Pi, gsl_vector* btwid, gsl_matrix* cov_mat, double nu)
{
	gsl_matrix_set_identity(Pi);
	double s = gsl_matrix_get(cov_mat,0,0);
	
	int nrows=cov_mat->size1;

	gsl_matrix_set(Pi,0,0,0.0);	
	
	double val;
	int i;
	for (i=1; i<nrows; i++)
	{
		val=(-1.0)*gsl_matrix_get(cov_mat,0,i)/s;
		gsl_matrix_set(Pi,i,0,val);
	}
		
	gsl_vector_set_zero(btwid);
	gsl_vector_set(btwid,0,1.0);

	double prefactor=sqrt(4*pi)*nu;

	for (i=1; i<nrows; i++)
	{
		val=prefactor*gsl_matrix_get(cov_mat,0,i)/s;
		gsl_vector_set(btwid,i,val);
	}
}
		
	
