#ifndef __PROJECTION_H_INCLUDED__
#define __PROJECTION_H_INCLUDED__

#include <gsl/gsl_matrix.h>

void gslvector_outer_product(gsl_vector* vec1, gsl_vector* vec2, gsl_matrix* prod);
void construct_projector(gsl_matrix* Pi, gsl_vector* btwid, gsl_matrix* cov_mat, double nu);

#endif // __PROJECTION_H_INCLUDED__
