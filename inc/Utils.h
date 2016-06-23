#ifndef __UTILS_H_INCLUDED__
#define __UTILS_H_INCLUDED__

#include <gsl/gsl_matrix.h>

const double pi=acos(-1.0);


void print_gslvec(gsl_vector* vec);
void print_gslmat(gsl_matrix* mat);
void linear_vector_ramp(gsl_vector * vec, double r_final);
void cholesky_decomp(gsl_matrix* M, gsl_matrix* A);
void print_eigenstuff(gsl_matrix* M);


#endif // __UTILS_H_INCLUDED__
