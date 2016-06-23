#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>

#include <stdexcept>
#include <iostream>

#include "Utils.h"
#include "CovGen.h"

double P(double k)
{
	//Using a decaying exponential as a dummy power spectrum.
	//FIXME: normalization?
	//return exp((-0.01)*k);
	return 1.0;
}


//Conditional sinc function with x=0 branch to handle removable
//discontinuity.
double ksinc(double r, double k, double epsilon)
{
	if(r<epsilon)
	{
		return k;
	}
	else
	{
		return sin(k*r)/r;
	}
}

//Integrand for covariance integral.
double f(double k, void* params)
{
	double r1=((double*)params)[0];
	double r2=((double*)params)[1];
	double epsilon=((double*)params)[2];
	return P(k)*ksinc(r1,k,epsilon)*ksinc(r2,k,epsilon);	
}			

double calculate_cov_element(double r1,double r2,double k_cutoff,double epsilon)
{

	std::cout << "R1: " << r1 << "  R2: " << r2 << "  epsilon: " << epsilon << std::endl;
	//Initialize this outside the function and pass pointer for speedup?
	gsl_integration_workspace *w
		= gsl_integration_workspace_alloc(10000);

	double params[3] = {r1,r2,epsilon};

	double result, error;
	gsl_function F;
	F.function = &f;
	F.params=params;
	
	gsl_integration_qag (&F, 0, k_cutoff, 0, 1e-6, 10000, 4, w, &result, &error);

	gsl_integration_workspace_free (w);

	return (16*pow(pi,2))*result;	
}


double precalc_calculate_cov_element(double r1,double r2,double k_cutoff,double epsilon)
{
	double result;

	if ((r1<epsilon)&&(r2<epsilon))
	{
		result = pow(k_cutoff,3)/3.0;
	}
	else if ((r1<epsilon)&&(r2>epsilon))
	{
		result = (1.0/pow(r2,2))*((sin(k_cutoff*r2)/r2) - k_cutoff*cos(k_cutoff*r2));
	}
	else if ((r1>epsilon)&&(r2<epsilon))
	{
		result = (1.0/pow(r1,2))*((sin(k_cutoff*r1)/r1) - k_cutoff*cos(k_cutoff*r1));
	}
	else
	{
		result = (1.0/(2*r1*r2))*(ksinc((r1-r2),k_cutoff,epsilon)-ksinc((r1+r2),k_cutoff,epsilon));
	}

	std::cout << "R1: " << r1 << "  R2: " << r2 << "  epsilon: " << epsilon << std::endl;
	std::cout << "RESULT:  " << result << std::endl;

	return (16*pow(pi,2))*result;	
}

void calculate_cov(gsl_matrix* cov_mat, gsl_vector* sample_radii, double k_cutoff)
{
	//Fill the covariance matrix by calculating the corresponding field
	//covariances with calculate_cov_element().


	//Check matrix and vector dimensions.
	if(cov_mat->size1 != cov_mat->size2)
	{
		throw std::runtime_error("Covariance matrix (gsl_matrix*) cov_mat is not square!");
	}
	if(cov_mat->size1 != sample_radii->size)
	{
		throw std::runtime_error("Covariance matrix (gsl_matrix*) cov_mat does not have the same dimension as (gsl_vector*) sample_radii!");
	}



	//Calculate an epsilon for use in zero-comparison for the r=0 branch
	//of the sinc function! This is a tad complicated because these 
	//calculations should work for nsamples >= 0.
	int nsamples=sample_radii->size;
	double epsilon;
	if (nsamples==1)
	{
		epsilon=1.0; //Case of origin sample only.
	}
	else
	{
		epsilon=gsl_vector_get(sample_radii,1)/2.0; //Multi-sample
	}

	int i,j;
	double r1,r2;
	double cov_element;
	for(i=0; i<nsamples; i++)
	{
		r1=gsl_vector_get(sample_radii,i);
		for(j=0; j<i+1; j++) //Calculate on bottom triangle only: symmetric matrix.
		{
			r2=gsl_vector_get(sample_radii,j);
			cov_element=calculate_cov_element(r1,r2,k_cutoff,epsilon);
			//cov_element=precalc_calculate_cov_element(r1,r2,k_cutoff,epsilon);
	
			gsl_matrix_set(cov_mat,i,j,cov_element);
		
			if (i!=j) //Symmetrize covariance matrix without overwriting the diagonal.
			{
				gsl_matrix_set(cov_mat,j,i,cov_element);
			}
		}
	}
}	

void precalculate_cov(gsl_matrix* cov_mat, gsl_vector* sample_radii, double k_cutoff)
{
	//Fill the covariance matrix by calculating the corresponding field
	//covariances with calculate_cov_element().


	//Check matrix and vector dimensions.
	if(cov_mat->size1 != cov_mat->size2)
	{
		throw std::runtime_error("Covariance matrix (gsl_matrix*) cov_mat is not square!");
	}
	if(cov_mat->size1 != sample_radii->size)
	{
		throw std::runtime_error("Covariance matrix (gsl_matrix*) cov_mat does not have the same dimension as (gsl_vector*) sample_radii!");
	}



	//Calculate an epsilon for use in zero-comparison for the r=0 branch
	//of the sinc function! This is a tad complicated because these 
	//calculations should work for nsamples >= 0.
	int nsamples=sample_radii->size;
	double epsilon;
	if (nsamples==1)
	{
		epsilon=1.0; //Case of origin sample only.
	}
	else
	{
		epsilon=gsl_vector_get(sample_radii,1)/2.0; //Multi-sample
	}

	int i,j;
	double r1,r2;
	double cov_element;
	for(i=0; i<nsamples; i++)
	{
		r1=gsl_vector_get(sample_radii,i);
		for(j=0; j<i+1; j++) //Calculate on bottom triangle only: symmetric matrix.
		{
			r2=gsl_vector_get(sample_radii,j);
			//cov_element=calculate_cov_element(r1,r2,k_cutoff,epsilon);
			cov_element=precalc_calculate_cov_element(r1,r2,k_cutoff,epsilon);
	
			gsl_matrix_set(cov_mat,i,j,cov_element);
		
			if (i!=j) //Symmetrize covariance matrix without overwriting the diagonal.
			{
				gsl_matrix_set(cov_mat,j,i,cov_element);
			}
		}
	}
}	
