#ifndef PYGMIN_POTENTIAL_H
#define PYGMIN_POTENTIAL_H

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

		virtual double get_energy(Array<double> x) { return 0.0; }
		virtual double get_energy_gradient(Array<double> x, Array<double> grad) { return 0.0; }

		/**
		 * compute the energy and gradient, but don't intialize the gradient to zero
		 */
		virtual double add_energy_gradient(Array<double> x, Array<double> grad) { return 0.0; }
	};
}

#endif
