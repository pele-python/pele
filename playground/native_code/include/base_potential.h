#ifndef PYGMIN_POTENTIAL_H
#define PYGMIN_POTENTIAL_H

#include <assert.h>
#include <vector>
#include <stdexcept>
#include <iostream>
#include "array.h"

namespace pele {
	/***
	 * basic potential interface for native potentials
	 */
	class BasePotential {
	public:
		virtual ~BasePotential() {}

		/**
		 * Return the energy of configuration x.  This is the only function which
		 * must be overloaded
		 */
		virtual double get_energy(Array<double> x) { return 0.0; }

		/**
		 * compute the energy and gradient, but don't intialize the gradient to zero
		 */
		virtual double add_energy_gradient(Array<double> x, Array<double> grad) { return 0.0; }

		/**
		 * compute the energy and gradient
		 */
		virtual double get_energy_gradient(Array<double> x, Array<double> grad) { return 0.; }

		/**
		 * compute the numerical gradient
		 */
		void numerical_gradient(Array<double> x, Array<double> grad, double eps=1e-6)
		{
			assert(x.size() == grad.size());
			for (size_t i=0; i<grad.size(); ++i){
				grad[i] = 0.;
			}

			Array<double> xnew(x.copy());
			for (size_t i=0; i<xnew.size(); ++i){
				xnew[i] -= eps;
				double eminus = get_energy(xnew);
				xnew[i] += 2. * eps;
				double eplus = get_energy(xnew);
				grad[i] = (eplus - eminus) / (2. * eps);
				xnew[i] = x[i];
			}
		}

	};
}

#endif
