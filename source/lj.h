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
#include "atomlist_potential.h"
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
        double inline energy(double r2, size_t atom_i, size_t atom_j) const {
            double ir2 = 1.0/r2;
            double ir6 = ir2*ir2*ir2;
            double ir12 = ir6*ir6;

            return -_C6*ir6 + _C12*ir12;
        }

        /* calculate energy and gradient from distance squared, gradient is in g/|rij| */
        double inline energy_gradient(double r2, double *gij, size_t atom_i, size_t atom_j) const {
            double ir2 = 1.0/r2;
            double ir6 = ir2*ir2*ir2;
            double ir12 = ir6*ir6;

            *gij = (_12C12 * ir12 - _6C6 * ir6) * ir2;
            return -_C6*ir6 + _C12*ir12;
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
     * Pairwise Lennard-Jones potential
     */
//    class LJAtomList : public AtomListPotential<lj_interaction, cartesian_distance> {
//        public:
//            LJAtomList(double C6, double C12, Array<size_t> atoms1, Array<size_t> atoms2) :
//                AtomListPotential<lj_interaction, cartesian_distance>(
//                        new lj_interaction(C6, C12),
//                        new cartesian_distance(),
//                        atoms1, atoms2)
//            {}
//    };

    class LJFrozen : public FrozenPotentialWrapper<LJ> {
        public:
            LJFrozen(double C6, double C12, Array<double> & reference_coords, Array<long int> & frozen_dof)
                : FrozenPotentialWrapper< LJ > 
                  ( new LJ(C6, C12), reference_coords, frozen_dof ) {}
    };

    /**
     * Pairwise Lennard-Jones potential with interaction lists
     */
    class LJNeighborList : public SimplePairwiseNeighborList< lj_interaction > {
        public:
            LJNeighborList(Array<long int> & ilist, double C6, double C12)
                :  SimplePairwiseNeighborList< lj_interaction > ( new lj_interaction(C6, C12), ilist) {}
    };
}
#endif
