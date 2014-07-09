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
#include <memory>

namespace pele {

/**
 * Pairwise interaction for lennard jones 
 */
struct lj_interaction {
    double const _C6, _C12;
    double const _6C6, _12C12;
    double const _42C6, _156C12;
    lj_interaction(double C6, double C12) 
        : _C6(C6), _C12(C12),
          _6C6(6.*_C6), _12C12(12.*_C12),
          _42C6(42.*_C6), _156C12(156.*_C12)
    {}

    /* calculate energy from distance squared */
    double inline energy(double r2, size_t atom_i, size_t atom_j) const 
    {
        double ir2 = 1.0/r2;
        double ir6 = ir2*ir2*ir2;
        double ir12 = ir6*ir6;

        return -_C6*ir6 + _C12*ir12;
    }

    /* calculate energy and gradient from distance squared, gradient is in g/|rij| */
    double inline energy_gradient(double r2, double *gij, size_t atom_i, size_t atom_j) const 
    {
        double ir2 = 1.0/r2;
        double ir6 = ir2*ir2*ir2;
        double ir12 = ir6*ir6;

        *gij = (_12C12 * ir12 - _6C6 * ir6) * ir2;
        return -_C6*ir6 + _C12*ir12;
    }

    double inline energy_gradient_hessian(double r2, double *gij, double *hij,
            size_t atom_i, size_t atom_j) const 
    {
        double ir2 = 1.0/r2;
        double ir6 = ir2*ir2*ir2;
        double ir12 = ir6*ir6;
        *gij = (_12C12 * ir12 - _6C6 * ir6) * ir2;
        *hij = (_156C12 * ir12 - _42C6 * ir6) * ir2;
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
            : SimplePairwisePotential< lj_interaction > (
                    std::make_shared<lj_interaction>(C6, C12) ) 
    {}
};

/**
 * Pairwise Lennard-Jones potential in a rectangular box
 */
class LJPeriodic : public SimplePairwisePotential< lj_interaction, periodic_distance<3> > {
public:
    LJPeriodic(double C6, double C12, Array<double> const boxvec)
        : SimplePairwisePotential< lj_interaction, periodic_distance<3>> (
                std::make_shared<lj_interaction>(C6, C12),
                std::make_shared<periodic_distance<3>>(boxvec)
                ) 
    {}
};

/**
 * Pairwise Lennard-Jones potential with frozen atoms
 */
class LJFrozen : public FrozenPotentialWrapper<LJ> {
public:
    LJFrozen(double C6, double C12, Array<double> & reference_coords, Array<size_t> & frozen_dof)
        : FrozenPotentialWrapper< LJ > (std::make_shared<LJ>(C6, C12),
                reference_coords, frozen_dof ) 
    {}
};

/**
 * Pairwise Lennard-Jones potential with interaction lists
 */
class LJNeighborList : public SimplePairwiseNeighborList< lj_interaction > {
public:
    LJNeighborList(Array<long int> & ilist, double C6, double C12)
        :  SimplePairwiseNeighborList<lj_interaction>(
                std::make_shared<lj_interaction>(C6, C12), ilist)
    {}
};
}
#endif
