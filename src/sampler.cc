#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>

#include <iostream>

using namespace std; //FIXME: use explicit namespaces in refined version.

const double pi=acos(-1.0);


void print_gslvec(gsl_vector* vec)
{
	int i;
	for(i=0; i<vec->size; i++)
	{
		cout << gsl_vector_get(vec,i) << endl;
	}
}

void print_gslmat(gsl_matrix* mat)
{
	cout <<"#--------------------------------------------------#"<<endl;
	int i,j;
	for(i=0; i<mat->size1; i++)
	{
		for(j=0; j<mat->size2; j++)
		{
			cout << gsl_matrix_get(mat,i,j) << "  ";
		}

		cout << endl;
	}
	cout <<"#--------------------------------------------------#"<<endl;
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

double P(double k)
{
	//Using a decaying exponential as a dummy power spectrum.
	//FIXME: normalization?
	return exp((-1.0)*k);
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

	//Initialize this outside the function and pass pointer for speedup?
	gsl_integration_workspace *w
		= gsl_integration_workspace_alloc(1000);

	double params[3] = {r1,r2,epsilon};

	double result, error;
	gsl_function F;
	F.function = &f;
	F.params=params;
	
	gsl_integration_qag (&F, 0, k_cutoff, 0, 1e-7, 1000, 4, w, &result, &error);

	gsl_integration_workspace_free (w);

	return (16*pow(pi,2))*result;	
}

void calculate_cov(gsl_matrix* cov_mat, gsl_vector* sample_radii, double k_cutoff)
{
	//Fill the covariance matrix by calculating the corresponding field
	//covariances with calculate_cov_element().

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
	
			gsl_matrix_set(cov_mat,i,j,cov_element);
		
			if (i!=j) //Symmetrize covariance matrix without overwriting the diagonal.
			{
				gsl_matrix_set(cov_mat,j,i,cov_element);
			}
		}
	}
}	

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
