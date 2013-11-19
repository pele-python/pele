#ifndef _PELE_WCA_H
#define _PELE_WCA_H

#include "simple_pairwise_potential.h"
#include "simple_pairwise_ilist.h"
#include "atomlist_potential.h"
#include "distance.h"

using std::exp;
using std::sqrt;

namespace pele {

    /**
     * Pairwise interaction for Weeks-Chandler-Andersen (WCA) potential
     */
    struct WCA_interaction {
        double const _C6, _C12;
        double const _6C6, _12C12;
        double const _coff, _eps; //cutoff distance for WCA potential

        WCA_interaction(double C6, double C12, double eps) :
            _C6(C6), _C12(C12),
            _6C6(6.*_C6), _12C12(12.*_C12),
        	_coff(pow(2.*C6,1./6)), _eps(eps)
        {}

        /* calculate energy from distance squared */
        double energy(double r2) const {
            double E;
        	double ir2 = 1.0/r2;
            double ir6 = ir2*ir2*ir2;
            double ir12 = ir6*ir6;
            if(sqrt(r2) < _coff)
            	E = _eps*(-_C6*ir6 + _C12*ir12 + 1.0/4);
            else
            	E = 0.;

            return E;
        }

        /* calculate energy and gradient from distance squared, gradient is in g/|rij| */
        double energy_gradient(double r2, double *gij) const {
        	double E;
			double ir2 = 1.0/r2;
			double ir6 = ir2*ir2*ir2;
			double ir12 = ir6*ir6;
			if(sqrt(r2) < _coff)
			{
				E = E = _eps*(-_C6*ir6 + _C12*ir12 + 1./4);
				*gij = _eps*(_12C12 * ir12 - _6C6 * ir6) * ir2;
			}
			else
			{
				E = 0.;
				*gij = 0;
			}

            return E;
        }
    };


    //
    // combine the components (interaction, looping method, distance function) into
    // defined classes
    //

    /**
     * Pairwise WCA potential
     */
    class WCA : public SimplePairwisePotential< WCA_interaction >
    {
        public:
            WCA(double C6, double C12, double eps)
                : SimplePairwisePotential< WCA_interaction > ( new WCA_interaction(C6, C12, eps) ) {}
    };

    /**
     * Pairwise WCA potential in a rectangular box
     */
    class WCAPeriodic : public SimplePairwisePotential< WCA_interaction, periodic_distance > {
        public:
            WCAPeriodic(double C6, double C12, double eps, double const *boxvec)
                : SimplePairwisePotential< WCA_interaction, periodic_distance> (
                        new WCA_interaction(C6, C12, eps),
                        new periodic_distance(boxvec[0], boxvec[1], boxvec[2])
                        )
            {}
    };

    /**
     * Pairwise WCA potential with interaction lists
     */
    class WCA_interaction_list : public SimplePairwiseInteractionList< WCA_interaction > {
        public:
            WCA_interaction_list(Array<long int> & ilist, double C6, double C12, double eps)
                :  SimplePairwiseInteractionList< WCA_interaction > ( new WCA_interaction(C6, C12, eps), ilist) {}
    };
}
#endif
