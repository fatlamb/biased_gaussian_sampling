#include "Spectra.h"
#include "Utils.h"

double gaussian(double x, double mu, double sigma, double scale)
{
	return scale*(1/(sigma*sqrt(2*pi)))*exp(-0.5*pow((x-mu)/sigma,2));
}
