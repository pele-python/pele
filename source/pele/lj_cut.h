#ifndef _PELE_LJ_CUT_H
#define _PELE_LJ_CUT_H

#include "simple_pairwise_potential.h"
#include "simple_pairwise_ilist.h"
#include "atomlist_potential.h"
#include "distance.h"
#include "frozen_atoms.h"

namespace pele {

/**
 * Define a pairwise interaction for lennard jones with a cutoff.  The
 * potential goes smoothly to zero using a second order
 * polynomial.
 */
struct lj_interaction_cut_smooth {
    double const _C6, _C12;
    double const _6C6, _12C12;
    double const _42C6, _156C12;
    double const _rcut2;
    double const _A0;
    double const _A2;
    double const _2A2;
    lj_interaction_cut_smooth(double C6, double C12, double rcut) 
        : 
            _C6(C6), _C12(C12), 
            _6C6(6.*_C6), _12C12(12.*_C12),
            _42C6(42.*_C6), _156C12(156.*_C12),
            _rcut2(rcut*rcut),
            //A0 = 4.0*(sig**6/rcut**6) - 7.0*(sig**12/rcut**12)
            _A0( 4.*_C6 / (_rcut2*_rcut2*_rcut2) - 7.*_C12/(_rcut2*_rcut2*_rcut2*_rcut2*_rcut2*_rcut2)),
            //A2 = (-3.0*(sig6/rcut**8) + 6.0*(sig12/rcut**14))
            _A2( -3.*_C6 / (_rcut2*_rcut2*_rcut2*_rcut2) + 6.*_C12/(_rcut2*_rcut2*_rcut2*_rcut2*_rcut2*_rcut2*_rcut2)),
            _2A2(2.*_A2)
    {}

    /* calculate energy from distance squared */
    double inline energy(double r2, size_t atom_i, size_t atom_j) const 
    {
        if (r2 >= _rcut2) {
            return 0.;
        }
        double ir2 = 1.0/r2;
        double ir6 = ir2*ir2*ir2;
        double ir12 = ir6*ir6;

        return -_C6*ir6 + _C12*ir12 + _A0 + _A2*r2;
    }

    /* calculate energy and gradient from distance squared, gradient is in g/|rij| */
    double inline energy_gradient(double r2, double *gij, size_t atom_i, size_t atom_j) const 
    {
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

    double inline energy_gradient_hessian(double r2, double *gij, double *hij, size_t atom_i, size_t atom_j) const 
    {
        if (r2 >= _rcut2) {
            *gij = 0.;
            *hij = 0.;
            return 0.;
        }
        double ir2 = 1.0/r2;
        double ir6 = ir2*ir2*ir2;
        double ir12 = ir6*ir6;

        *gij = (_12C12 * ir12 - _6C6 * ir6) * ir2 - _2A2;
        *hij = (_156C12 * ir12 - _42C6 * ir6) * ir2 + _2A2;
        return -_C6*ir6 + _C12*ir12 + _A0 + _A2*r2;
    }
};

/**
 * Pairwise Lennard-Jones potential with smooth cutoff
 */
class LJCut : public SimplePairwisePotential< lj_interaction_cut_smooth > {
public:
    LJCut(double C6, double C12, double rcut)
        : SimplePairwisePotential< lj_interaction_cut_smooth > (
                std::make_shared<lj_interaction_cut_smooth>(C6, C12, rcut) ) 
    {}
};

/**
 * Pairwise Lennard-Jones potential with smooth cutoff in a rectangular box
 */
class LJCutPeriodic : public SimplePairwisePotential< lj_interaction_cut_smooth, periodic_distance<3>> {
public:
    LJCutPeriodic(double C6, double C12, double rcut, Array<double> const boxvec)
        : SimplePairwisePotential< lj_interaction_cut_smooth, periodic_distance<3>> (
                std::make_shared<lj_interaction_cut_smooth>(C6, C12, rcut),
                std::make_shared<periodic_distance<3>>(boxvec)
                ) 
    {}
};

/**
 * Pairwise Lennard-Jones potential with smooth cutoff with loops done
 * using atom lists
 */
class LJCutAtomList : public AtomListPotential<lj_interaction_cut_smooth, cartesian_distance<3>> {
public:
LJCutAtomList(double C6, double C12, double rcut, Array<size_t> atoms1, Array<size_t> atoms2) 
    : AtomListPotential<lj_interaction_cut_smooth, cartesian_distance<3>>(
            std::make_shared<lj_interaction_cut_smooth>(C6, C12, rcut),
            std::make_shared<cartesian_distance<3>>(), atoms1, atoms2)
{}

LJCutAtomList(double C6, double C12, double rcut, Array<size_t> atoms1) 
    : AtomListPotential<lj_interaction_cut_smooth, cartesian_distance<3>>(
            std::make_shared<lj_interaction_cut_smooth>(C6, C12, rcut),
            std::make_shared<cartesian_distance<3>>(), atoms1)
{}
};

/**
 * Pairwise Lennard-Jones potential with smooth cutoff with loops done
 * using atom lists
 */
class LJCutPeriodicAtomList : public AtomListPotential<lj_interaction_cut_smooth, periodic_distance<3>> {
public:
    LJCutPeriodicAtomList(double C6, double C12, double rcut, pele::Array<double> boxvec,
            Array<size_t> atoms1, Array<size_t> atoms2) 
        : AtomListPotential<lj_interaction_cut_smooth, periodic_distance<3> >(
                std::make_shared<lj_interaction_cut_smooth>(C6, C12, rcut),
                std::make_shared<periodic_distance<3>>(boxvec), atoms1,
                atoms2)
    {}

    LJCutPeriodicAtomList(double C6, double C12, double rcut, pele::Array<double> boxvec,
            Array<size_t> atoms1) 
        : AtomListPotential<lj_interaction_cut_smooth, periodic_distance<3>>(
                std::make_shared<lj_interaction_cut_smooth>(C6, C12, rcut),
                std::make_shared<periodic_distance<3>>(boxvec), atoms1)
    {}
};

}

#endif
