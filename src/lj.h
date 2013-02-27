/******************************+
 * This is an example on how to define a new c potential. The class
 * internally calls the lennard jones fortran function but manages all
 * the parameters.
 *
 * The python bindings for this potential are in lj.pyx
 *
 * For an alternative pure cython implementation for this interface
 * see _lj_cython.pyx
 *
 */

#ifndef PYGMIN_LJ_H
#define PYGMIN_LJ_H

// subroutine ljenergy_gradient( coords, natoms, e, grad, eps, sig, periodic, boxl )


#include "potential.h"

namespace pygmin {

	class LJ  : public Potential
	{
		double _eps, _sigma;
	public:
		LJ(double sigma, double eps) : _eps(eps), _sigma(sigma) {}

		double get_energy(Array &x);
		double get_energy_gradient(Array &x, Array &out);
	};

}

#endif
