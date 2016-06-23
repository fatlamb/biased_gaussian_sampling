#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <math.h>

#include <stdexcept>
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
    //Check matrix and vector dimensions.
    if(cov_mat->size1 != cov_mat->size2)
    {
        throw std::runtime_error("Covariance matrix (gsl_matrix*) cov_mat is not square!");
    }
    if(Pi->size1 != Pi->size2)
    {
        throw std::runtime_error("Projection matrix (gsl_matrix*) Pi is not square!");
    }

    if(cov_mat->size1 != Pi->size1)
    {
        throw std::runtime_error("Covariance matrix (gsl_matrix*) cov_mat does not have the same dimensions as projection matrix (gsl_matrix*) Pi");
    }
    if(cov_mat->size1 != btwid->size)
    {
        throw std::runtime_error("Covariance matrix (gsl_matrix*) cov_mat does not have the same dimension as offset vector (gsl_vector*) btwid");
    }


	gsl_matrix_set_identity(Pi);
	double s = gsl_matrix_get(cov_mat,0,0);
	
	int n=cov_mat->size1;

	gsl_matrix_set(Pi,0,0,0.0);	
	
	double val;
	int i;
	for (i=1; i<n; i++)
	{
		val=(-1.0)*gsl_matrix_get(cov_mat,0,i)/s;
		gsl_matrix_set(Pi,i,0,val);
	}
		
	gsl_vector_set_zero(btwid);

	double prefactor=sqrt(4*pi)*nu;

	gsl_vector_set(btwid,0,prefactor);


	for (i=1; i<n; i++)
	{
		val=prefactor*gsl_matrix_get(cov_mat,0,i)/s;
		gsl_vector_set(btwid,i,val);
	}
}
		
	
