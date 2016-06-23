#include <gsl/gsl_matrix.h>

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "WriteOut.h"

void write_gslvecs(std::vector<gsl_vector*> vecs, std::string filename)
{
    //Check vector dimensions.
	bool mismatch=false;

	int i;
	for (i=1; i<vecs.size(); i++)
	{
		if(vecs[i]->size != vecs[0]->size)
		{
			mismatch=true;
			break;
		}
	}

    if(mismatch)
    {
        throw std::runtime_error("Mismatched (gsl_vector*) dimensions in std::vector<gsl_vector*> vecs");
    }


	std::ofstream file;
	file.open(filename.c_str());
	
	int j,n,nvecs;
	n=vecs[0]->size;
	nvecs=vecs.size();

	std::stringstream ss (std::stringstream::in | std::stringstream::out);
	std::string word;

	for (i=0; i<n; i++)
	{
		for(j=0; j<nvecs; j++)
		{
			ss << gsl_vector_get(vecs[j],i);
			word = ss.str();

			ss.clear();//clear any bits set
			ss.str(std::string());

			file << word;

			if (j==(nvecs-1))
			{
				file << "\n";
			}
			else
			{
				file << " ";
			}
		}
	}

	file.close();	
}
