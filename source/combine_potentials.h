#ifndef _PELE_COMBINE_POTENTIALS_H_
#define _PELE_COMBINE_POTENTIALS_H_
#include <list>
#include "base_potential.h"

namespace pele {

	/**
	 * Potential wrapper which wraps multiple potentials to
	 * act as one potential.  This can be used to implement
	 * multiple types of interactions, e.g. a system with two
	 * types atoms
	 */
	class CombinedPotential : public BasePotential{
	protected:
		std::list<BasePotential*> _potentials;

	public:
		CombinedPotential() {}

		/**
		 * destructor: destroy all the potentials in the list
		 */
		virtual ~CombinedPotential()
		{
			BasePotential * potential_ptr;
			std::list<BasePotential*>::iterator iter;
			for (iter = _potentials.begin(); iter != _potentials.end(); ++iter){
				potential_ptr = *iter;
				if (potential_ptr != NULL) {
					delete potential_ptr;
				}
			}
		}

		/**
		 * add a potential to the list
		 */
		virtual void add_potential(BasePotential * potential)
		{
			_potentials.push_back(potential);
		}

		virtual double get_energy(Array<double> x)
		{
			double energy = 0.;
			std::list<BasePotential*>::iterator iter;
			for (iter = _potentials.begin();
					iter != _potentials.end(); ++iter){
				energy += (*iter)->get_energy(x);
			}
			return energy;
		}

		virtual double get_energy_gradient(Array<double> x, Array<double> grad)
		{
			double energy = 0.;
			for (size_t i=0; i<grad.size(); ++i){
				grad[i] = 0;
			}

			std::list<BasePotential*>::iterator iter;
			for (iter = _potentials.begin();
					iter != _potentials.end(); ++iter){
				energy += (*iter)->add_energy_gradient(x, grad);
			}
			return energy;
		}

	};
}

#endif
