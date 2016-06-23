#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

#include <stdexcept>
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

void cholesky_decomp(gsl_matrix* M, gsl_matrix* A)
{

    //Check matrix dimensions.
    if(M->size1 != M->size2)
    {
        throw std::runtime_error("(gsl_matrix*) M is not square!");
    }
    if(A->size1 != A->size2)
    {
        throw std::runtime_error("(gsl_matrix*) A is not square!");
    }
    if(M->size1 != A->size1)
    {
        throw std::runtime_error("(gsl_matrix*) M dimensions do not match those of (gsl_matrix*) A!");
    }



	int ret; //What do do with ret? Error code?
	ret = gsl_linalg_cholesky_decomp(M);

	gsl_matrix_set_zero(A);


	int n=M->size1;

	int i,j;
	for (i=0; i<n; i++)
	{
		for (j=0; j<i+1; j++)
		{
			gsl_matrix_set(A,i,j,gsl_matrix_get(M,i,j));
		}
	}			

}


void print_eigenstuff(gsl_matrix* M)
{
	
	//Check matrix dimensions.
	if(M->size1 != M->size2)
	{
		throw std::runtime_error("(gsl_matrix*) M is not square!");
 	}
	
	int n = M->size1;

	gsl_vector *eval = gsl_vector_alloc (n);
	gsl_matrix *evec = gsl_matrix_alloc (n, n);
	
	gsl_eigen_symmv_workspace * w = 
		gsl_eigen_symmv_alloc (n);
	
	gsl_eigen_symmv (M, eval, evec, w);
	
	gsl_eigen_symmv_free (w);
	
	gsl_eigen_symmv_sort (eval, evec, 
		GSL_EIGEN_SORT_ABS_ASC);
	
	int i;
	for (i = 0; i < n; i++)
	{
		double eval_i 
		 	= gsl_vector_get (eval, i);
		gsl_vector_view evec_i 
			= gsl_matrix_column (evec, i);
		
		printf ("eigenvalue = %g\n", eval_i);
		//printf ("eigenvector = \n");
		//gsl_vector_fprintf (stdout, 
		//&evec_i.vector, "%g");
	}
	
	gsl_vector_free (eval);
	gsl_matrix_free (evec);
}
