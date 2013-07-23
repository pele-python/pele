#include <iostream>
#include "lj.h"

extern "C" {
// this is in fortran
extern void ljenergy_gradient_(double *x, int *natoms, double *e, double *grad,
		double *eps, double *sig, int *peridic, double *boxl);
};

namespace pele {

double LJ::get_energy(Array &x)
{
	return 0.;
}

double LJ::get_energy_gradient(Array &x, Array &grad)
{
	int natoms = x.size()/3;
	double e;
	int periodic = 0;
	double boxl = -1;

	ljenergy_gradient_(x.data(), &natoms, &e, grad.data(),
			&_eps, &_sigma, &periodic, &boxl);

	return e;
}

}
