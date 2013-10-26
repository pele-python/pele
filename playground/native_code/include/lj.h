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

#include "simple_pairwise_potential.h"
#include "simple_pairwise_ilist.h"
#include "distance.h"
#include "frozen_atoms.h"

namespace pele {

    /**
   * Pairwise interaction for lennard jones 
   */
    struct lj_interaction {
        double const _C6, _C12;
        double const _6C6, _12C12;
        lj_interaction(double C6, double C12) : 
            _C6(C6), _C12(C12),
            _6C6(6.*_C6), _12C12(12.*_C12) 
        {}

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

            *gij = (_12C12 * ir12 - _6C6 * ir6) * ir2;
            return -_C6*ir6 + _C12*ir12;
        }
    };

    /**
   * Define a pairwise interaction for lennard jones with a cutoff.  The
   * potential goes smoothly to zero using a second order
   * polynomial.
   */
    struct lj_interaction_cut_smooth {
        double const _C6, _C12;
        double const _6C6, _12C12;
        double const _rcut2;
        double const _A0;
        double const _A2;
        double const _2A2;
        lj_interaction_cut_smooth(double C6, double C12, double rcut) 
            : 
                _C6(C6), _C12(C12), 
                _6C6(6.*_C6), _12C12(12.*_C12),
                _rcut2(rcut*rcut),
                //A0 = 4.0*(sig**6/rcut**6) - 7.0*(sig**12/rcut**12)
                _A0( 4.*_C6 / (_rcut2*_rcut2*_rcut2) - 7.*_C12/(_rcut2*_rcut2*_rcut2*_rcut2*_rcut2*_rcut2)),
                //A2 = (-3.0*(sig6/rcut**8) + 6.0*(sig12/rcut**14))
                _A2( -3.*_C6 / (_rcut2*_rcut2*_rcut2*_rcut2) + 6.*_C12/(_rcut2*_rcut2*_rcut2*_rcut2*_rcut2*_rcut2*_rcut2)),
                _2A2(2.*_A2)
        {}

        /* calculate energy from distance squared */
        double energy(double r2) const {
            if (r2 >= _rcut2) {
                return 0.;
            }
            double ir2 = 1.0/r2;
            double ir6 = ir2*ir2*ir2;
            double ir12 = ir6*ir6;

            return -_C6*ir6 + _C12*ir12 + _A0 + _A2*r2;
        }

        /* calculate energy and gradient from distance squared, gradient is in g/|rij| */
        double energy_gradient(double r2, double *gij) const {
            if (r2 >= _rcut2) {
                *gij = 0.;
                return 0.;
            }
            double ir2 = 1.0/r2;
            double ir6 = ir2*ir2*ir2;
            double ir12 = ir6*ir6;

            *gij = (_12C12 * ir12 - _6C6 * ir6) * ir2 - _2A2;
            return -_C6*ir6 + _C12*ir12 + _A0 + _A2*r2;
        }
    };

    //
    // combine the components (interaction, looping method, distance function) into
    // defined classes
    //

    /**
     * Pairwise Lennard-Jones potential
     */
    class LJ : public SimplePairwisePotential< lj_interaction > {
        public:
            LJ(double C6, double C12)
                : SimplePairwisePotential< lj_interaction > ( new lj_interaction(C6, C12) ) {}
    };

    /**
     * Pairwise Lennard-Jones potential in a rectangular box
     */
    class LJPeriodic : public SimplePairwisePotential< lj_interaction, periodic_distance > {
        public:
            LJPeriodic(double C6, double C12, double const *boxvec)
                : SimplePairwisePotential< lj_interaction, periodic_distance> ( 
                        new lj_interaction(C6, C12), 
                        new periodic_distance(boxvec[0], boxvec[1], boxvec[2])
                        ) 
            {}
    };

    /**
     * Pairwise Lennard-Jones potential with smooth cutoff
     */
    class LJCut : public SimplePairwisePotential< lj_interaction_cut_smooth > {
        public:
            LJCut(double C6, double C12, double rcut)
                : SimplePairwisePotential< lj_interaction_cut_smooth > ( new  lj_interaction_cut_smooth(C6, C12, rcut) ) {}
    };

    /**
     * Pairwise Lennard-Jones potential with smooth cutoff in a rectangular box
     */
    class LJCutPeriodic : public SimplePairwisePotential< lj_interaction_cut_smooth, periodic_distance > {
        public:
            LJCutPeriodic(double C6, double C12, double rcut, double const *boxvec)
                : SimplePairwisePotential< lj_interaction_cut_smooth, periodic_distance> ( 
                        new lj_interaction_cut_smooth(C6, C12, rcut), 
                        new periodic_distance(boxvec[0], boxvec[1], boxvec[2])
                        ) 
            {}
    };

    class LJFrozen : public FrozenPotentialWrapper<LJ> {
        public:
            LJFrozen(double C6, double C12, Array<double> & reference_coords, Array<long int> & frozen_dof)
                : FrozenPotentialWrapper< LJ > 
                  ( new LJ(C6, C12), reference_coords, frozen_dof ) {}
    };


    /**
     * Pairwise Lennard-Jones potential with interaction lists
     */
    class LJ_interaction_list : public SimplePairwiseInteractionList< lj_interaction > {
        public:
            LJ_interaction_list(Array<long int> & ilist, double C6, double C12)
                :  SimplePairwiseInteractionList< lj_interaction > ( new lj_interaction(C6, C12), ilist) {}
    };
}
#endif
