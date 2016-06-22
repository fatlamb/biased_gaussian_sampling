#ifndef __UTILS_H_INCLUDED__
#define __UTILS_H_INCLUDED__

#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>

const double pi=acos(-1.0);


void print_gslvec(gsl_vector* vec);
void print_gslmat(gsl_matrix* mat);
void linear_vector_ramp(gsl_vector * vec, double r_final);


#endif // __UTILS_H_INCLUDED__
