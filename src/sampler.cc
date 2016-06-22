#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include <iostream>

using namespace std; //FIXME: use explicit namespaces in refined version.

double pi=acos(-1.0);


double P(double k)
{
	//Using a decaying exponential as a dummy power spectrum.
	//FIXME: normalization?
	return exp((-1.0)*k);
}


//Integrand for covariance integral.
double f(double k, void* params)
{
	double r1=((double*)params)[0];
	double r2=((double*)params)[1];
	return P(k)*sin(r1*k)*sin(r2*k);	
}			

double calculate_cov(double r1,double r2,double k_cutoff)
{

	//Initialize this outside the function and pass pointer for speedup?
	gsl_integration_workspace *w
		= gsl_integration_workspace_alloc(1000);

	double radii[2] = {r1,r2};

	double result, error;
	gsl_function F;
	F.function = &f;
	F.params=radii;
	
	gsl_integration_qag (&F, 0, k_cutoff, 0, 1e-7, 1000, 4, w, &result, &error);



	printf ("result          = % .18f\n", result);
	//printf ("exact result    = % .18f\n", expected);
	printf ("estimated error = % .18f\n", error);
	//printf ("actual error    = % .18f\n", result - expected);
	printf ("intervals       = %zu\n", w->size);

	gsl_integration_workspace_free (w);

	return result;	
}

int main()
{
	cout << calculate_cov(1,2,10) << endl;
	return 0;
}
