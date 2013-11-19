#ifndef _PELE_HS_WCA_H
#define _PELE_HS_WCA_H

#include "simple_pairwise_potential.h"
#include "simple_pairwise_ilist.h"
#include "atomlist_potential.h"
#include "distance.h"

using std::exp;
using std::sqrt;

namespace pele {

    /**
     * Pairwise interaction for Hard Sphere + Weeks-Chandler-Andersen (HS_WCA) potential
     */
    struct HS_WCA_interaction {
        double const _C6, _C12;
        double const _6C6, _12C12;
        double const _coff, _infty; //cutoff distance for HS_WCA potential
        double const _eps;

        HS_WCA_interaction(double C6, double C12, double eps) :
            _C6(C6), _C12(C12),
            _6C6(6.*_C6), _12C12(12.*_C12),
        	_coff(pow(2.*C6,1./6)),
        	_infty(pow(10.0,20)), _eps(eps)
        {}

        /* calculate energy from distance squared */
        double energy(double r2, double r0) const {
            double E;
        	double r = sqrt(r2);
            double ir2 = 1.0/r2;
            double ir6 = ir2*ir2*ir2;
            double ir12 = ir6*ir6;
            if (r < r0)
            	E = _infty;
            else if(r < _coff)
            	E = _eps*(-_C6*ir6 + _C12*ir12 + 1.0/4);
            else
            	E = 0.;

            return E;
        }

        /* calculate energy and gradient from distance squared, gradient is in g/|rij| */
        double energy_gradient(double r2, double r0, double *gij) const {
        	double E;
        	double r = sqrt(r2);
			double ir2 = 1.0/r2;
			double ir6 = ir2*ir2*ir2;
			double ir12 = ir6*ir6;
			if (r < r0)
			{
				E = _infty;
				*gij = _infty;
			}
			else if(r < _coff)
			{
				E = _eps*(-_C6*ir6 + _C12*ir12 + 1.0/4);
				*gij = _eps*(_12C12 * ir12 - _6C6 * ir6) * ir2;
			}
			else
			{
				E = 0.;
				*gij = 0.;
			}

            return E;
        }
    };


    //
    // combine the components (interaction, looping method, distance function) into
    // defined classes
    //

    /**
     * Pairwise HS_WCA potential
     */
    class HS_WCA : public SimplePairwisePotential< HS_WCA_interaction >
    {
        public:
            HS_WCA(double C6, double C12, double eps)
                : SimplePairwisePotential< HS_WCA_interaction > ( new HS_WCA_interaction(C6, C12, eps) ) {}
    };

    /**
     * Pairwise HS_WCA potential in a rectangular box
     */
    class HS_WCAPeriodic : public SimplePairwisePotential< HS_WCA_interaction, periodic_distance > {
        public:
            HS_WCAPeriodic(double C6, double C12, double eps, double const *boxvec)
                : SimplePairwisePotential< HS_WCA_interaction, periodic_distance> (
                        new HS_WCA_interaction(C6, C12, eps),
                        new periodic_distance(boxvec[0], boxvec[1], boxvec[2])
                        )
            {}
    };

    /**
     * Pairwise HS_WCA potential with interaction lists
     */
    class HS_WCA_interaction_list : public SimplePairwiseInteractionList< HS_WCA_interaction > {
        public:
            HS_WCA_interaction_list(Array<long int> & ilist, double C6, double C12, double eps)
                :  SimplePairwiseInteractionList< HS_WCA_interaction > ( new HS_WCA_interaction(C6, C12, eps), ilist) {}
    };
}
#endif
