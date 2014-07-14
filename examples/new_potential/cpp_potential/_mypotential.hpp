/******************************+
 * This is an example on how to define a new c++ potential. 
 *
 * The python bindings for this potential are in mypotential.pyx
 *
 */


#include <memory>
#include "pele/simple_pairwise_potential.h"
#include "pele/distance.h"

namespace pele {

    /**
     * Pairwise interaction for modified lennard jones 
     */
    struct mypot_interaction {
        double const _sig12;
        double const _sig24;
        double const _eps;
        mypot_interaction(double sig, double eps) : 
            _sig12(pow(sig, 12)), 
            _sig24(_sig12 * _sig12),
            _eps(eps)
        {}

        /* calculate energy from distance squared */
        double inline energy(double r2, size_t atom_i, size_t atom_j) const {
            double ir2 = 1.0/r2;
            double ir6 = ir2*ir2*ir2;
            double ir12 = ir6*ir6;
            double ir24 = ir12*ir12;

            return 4. * _eps * (_sig24 * ir24 - _sig12 * ir12);
        }

        /* calculate energy and gradient from distance squared, gradient is in |g|/|rij| */
        double inline energy_gradient(double r2, double *gij, size_t atom_i, size_t atom_j) const {
            double ir2 = 1.0/r2;
            double ir6 = ir2*ir2*ir2;
            double ir12 = ir6*ir6;
            double ir24 = ir12*ir12;

            *gij = 4. * _eps * (24. * _sig24 * ir24 - 12. * _sig12 * ir12) * ir2;
            double energy = 4. * _eps * (_sig24 * ir24 - _sig12 * ir12);
            return energy;
        }

        double inline energy_gradient_hessian(double r2, double *gij, double *hij, size_t atom_i, size_t atom_j) const {
            double ir2 = 1.0/r2;
            double ir6 = ir2*ir2*ir2;
            double ir12 = ir6*ir6;
            double ir24 = ir12*ir12;

            double energy = 4. * _eps * (_sig24 * ir24 - _sig12 * ir12);
            *gij = 4. * _eps * (24. * _sig24 * ir24 - 12. * _sig12 * ir12) * ir2;
            *hij = 4. * _eps * (24. * 25. * _sig24 * ir24 - 12. * 13 * _sig12 * ir12) * ir2;
            return energy;
        }
    };


    //
    //

    /**
     * We now combine the components (interaction, looping method, distance
     * function) into defined classes which is the c++ potential class
     */
    class MyPot : public pele::SimplePairwisePotential<mypot_interaction, cartesian_distance<3> > {
        public:
            MyPot(double sig, double eps)
                : SimplePairwisePotential<mypot_interaction, cartesian_distance<3> > (std::make_shared<mypot_interaction>(sig, eps)) {}
    };
}
