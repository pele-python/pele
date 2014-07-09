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
    double const _6C6, _12C12, _42C6, _156C12;
    double const _coff2, _eps; //cutoff distance for WCA potential

    WCA_interaction(double sig, double eps) 
        : _C6(sig*sig*sig*sig*sig*sig),
          _C12(_C6*_C6), _6C6(6.*_C6),
          _12C12(12.*_C12), _42C6(42*_C6),
          _156C12(156*_C12),
          _coff2(pow(2.,1./3)*sig*sig), _eps(eps)
    {}

    /* calculate energy from distance squared */
    double energy(double r2, size_t atom_i, size_t atom_j) const 
    {
        double E;
        double ir2 = 1.0/r2;
        double ir6 = ir2*ir2*ir2;
        double ir12 = ir6*ir6;
        if (r2 < _coff2)
            E = 4.*_eps*(-_C6*ir6 + _C12*ir12) + _eps;
        else
            E = 0.;

        return E;
    }

    /* calculate energy and gradient from distance squared, gradient is in -(dv/drij)/|rij| */
    double energy_gradient(double r2, double *gij, size_t atom_i, size_t atom_j) const 
    {
        double E;
        double ir2 = 1.0/r2;
        double ir6 = ir2*ir2*ir2;
        double ir12 = ir6*ir6;
        if (r2 < _coff2) {
            E = 4.*_eps*(-_C6*ir6 + _C12*ir12) + _eps;
            *gij = 4.*_eps*(- _6C6 * ir6 + _12C12 * ir12) * ir2;
        } else {
            E = 0.;
            *gij = 0;
        }

        return E;
    }

    double inline energy_gradient_hessian(double r2, double *gij, double *hij, size_t atom_i, size_t atom_j) const 
    {
        double E;
        double ir2 = 1.0/r2;
        double ir6 = ir2*ir2*ir2;
        double ir12 = ir6*ir6;

        if (r2 < _coff2) {
            E = 4.*_eps*(-_C6*ir6 + _C12*ir12) + _eps;
            *gij = 4.*_eps*(- _6C6 * ir6 + _12C12 * ir12) * ir2;
            *hij = 4.*_eps*(- _42C6 * ir6 + _156C12 * ir12) * ir2;
        } else {
            E = 0.;
            *gij = 0;
            *hij=0;
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
class WCA : public SimplePairwisePotential< WCA_interaction, cartesian_distance<3>> {
public:
    WCA(double sig, double eps)
        : SimplePairwisePotential< WCA_interaction, cartesian_distance<3>>(
                std::make_shared<WCA_interaction>(sig, eps),
                std::make_shared<cartesian_distance<3>>()
                )
    {}
};

class WCA2D : public SimplePairwisePotential< WCA_interaction, cartesian_distance<2> > {
public:
    WCA2D(double sig, double eps)
        : SimplePairwisePotential< WCA_interaction, cartesian_distance<2>> (
                std::make_shared<WCA_interaction>(sig, eps),
                std::make_shared<cartesian_distance<2>>()
        ) 
    {}
};

/**
 * Pairwise WCA potential in a rectangular box
 */
class WCAPeriodic : public SimplePairwisePotential< WCA_interaction, periodic_distance<3> > {
public:
    WCAPeriodic(double sig, double eps, Array<double> const boxvec)
        : SimplePairwisePotential< WCA_interaction, periodic_distance<3>> (
                std::make_shared<WCA_interaction>(sig, eps),
                std::make_shared<periodic_distance<3>>(boxvec)
                )
    {}
};

/**
 * Pairwise WCA potential in a rectangular box
 */
class WCAPeriodic2D : public SimplePairwisePotential< WCA_interaction, periodic_distance<2> > {
public:
    WCAPeriodic2D(double sig, double eps, Array<double> const boxvec)
        : SimplePairwisePotential< WCA_interaction, periodic_distance<2>> (
                std::make_shared<WCA_interaction>(sig, eps),
                std::make_shared<periodic_distance<2>>(boxvec)
                )
    {}
};

/**
 * Pairwise WCA potential with interaction lists
 */
class WCANeighborList : public SimplePairwiseNeighborList< WCA_interaction > {
public:
    WCANeighborList(Array<long int> & ilist, double sig, double eps)
        :  SimplePairwiseNeighborList< WCA_interaction > (
                std::make_shared<WCA_interaction>(sig, eps), ilist)
    {}
};
}
#endif


