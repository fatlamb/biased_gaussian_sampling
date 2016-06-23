#ifndef __WRITEOUT_H_INCLUDED__
#define __WRITEOUT_H_INCLUDED__

#include <gsl/gsl_matrix.h>
#include <string>
#include <vector>

void write_gslvecs(std::vector<gsl_vector*> vecs, std::string filename);

#endif // __WRITEOUT_H_INCLUDED__
