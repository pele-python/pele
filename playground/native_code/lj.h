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

#include "include/simple_pairwise_potential.h"
#include "include/simple_pairwise_ilist.h"
namespace pele {

	/* define a pairwise interaction for lennard jones */
	struct lj_interaction {
		double _C6, _C12;
		lj_interaction(double C6, double C12) : _C6(C6), _C12(C12) {}

		/* calculate energy from distance squared */
		double energy(double r2) const {
			double ir2 = 1.0/r2;
			double ir6 = ir2*ir2*ir2;
			double ir12 = ir6*ir6;

			return -_C6*ir6 + _C12*ir12;
		}

		/* calculate energy and gradient from distance squared, gradient is in g/|rij| */
		double energy_gradient(double r2, double *gij) const {
			double ir2 = 1.0/r2;
			double ir6 = ir2*ir2*ir2;
			double ir12 = ir6*ir6;

			*gij = (12.0 * _C12 * ir12 -  6.0 * _C6 * ir6) * ir2;
			return -_C6*ir6 + _C12*ir12;
		}

	};

	// define lennard jones potential as a pairwise interaction
	class LJ : public SimplePairwisePotential< lj_interaction > {
	public:
		LJ(double C6, double C12)
			: SimplePairwisePotential< lj_interaction > ( new lj_interaction(C6, C12) ) {}
	};

	// define lennard jones potential as a pairwise interaction
	class LJ_interaction_list : public SimplePairwiseInteractionList< lj_interaction > {
	public:
		LJ_interaction_list(Array<int> & ilist, double C6, double C12)
			:  SimplePairwiseInteractionList< lj_interaction > ( new lj_interaction(C6, C12), ilist) {}
	};
}
#endif
