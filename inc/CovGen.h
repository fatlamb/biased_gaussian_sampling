#ifndef __COVGEN_H_INCLUDED__
#define __COVGEN_H_INCLUDED__

#include <gsl/gsl_matrix.h>

double P(double k);
double ksinc(double r, double k, double epsilon);
double f(double k, void* params);
double calculate_cov_element(double r1,double r2,double k_cutoff,double epsilon);
void calculate_cov(gsl_matrix* cov_mat, gsl_vector* sample_radii, double k_cutoff);
void precalculate_cov(gsl_matrix* cov_mat, gsl_vector* sample_radii, double k_cutoff);


#endif // __COVGEN_H_INCLUDED__
