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
	class Potential {
	public:
		virtual ~Potential() {}

		virtual double get_energy(Array x) {} ;
		virtual double get_energy_gradient(Array x, Array grad) {} ;
	};
}

#endif
