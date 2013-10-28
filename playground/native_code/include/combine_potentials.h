#ifndef _PELE_COMBINE_POTENTIALS_H_
#define _PELE_COMBINE_POTENTIALS_H_
#include <list>
#include "base_potential.h"

namespace pele {
	class CombinedPotential : public BasePotential{
	protected:
		std::list<BasePotential*> _potentials;

	public:
		CombinedPotential() {}
//		~CombinedPotential() {
//			// deallocate all the potentials
//			for (std::list<BasePotential>::iterator iter = _potentials.begin();
//					iter != _potentials.end(); ++iter){
//				delete *iter;
//			}
//
//		}


		void add_potential(BasePotential * potential)
		{
			_potentials.push_back(potential);
		}

		double get_energy(Array<double> x)
		{
			double energy = 0.;
			std::list<BasePotential*>::iterator iter;
			for (iter = _potentials.begin();
					iter != _potentials.end(); ++iter){
				energy += (*iter)->get_energy(x);
			}
			return energy;
		}

		double get_energy_gradient(Array<double> x, Array<double> grad)
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
