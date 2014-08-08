#ifndef _PELE_INVERSEPOWER_H
#define _PELE_INVERSEPOWER_H

#include "simple_pairwise_potential.h"
#include "simple_pairwise_ilist.h"
#include "atomlist_potential.h"
#include "distance.h"
#include <memory>

namespace pele {

/**
 * Pairwise interaction for Inverse power potential eps/pow * (1 - r/r0)^pow
 * the radii array allows for polydispersity
 * The most common exponents are:
 * pow=2   -> Hookean interaction
 * pow=2.5 -> Hertzian interaction
 *
 * Comments about this implementation:
 *
 * This implementation is using STL pow for all exponents. Performance wise this may not be
 * the fastest possible implementation, for instance using exp(pow*log) could be faster
 * (twice as fast according to some blogs). Another possibility is to make the power a template
 * parameter and then write template specialisations for the most common exponents, this however
 * makes the interface with python ugly because cython does not deal well with template
 * integer/double parameter yet.
 *
 * So in the interest of simplicity I would keep the implementation as is for now, maybe we
 * should consider a function that calls exp(pow*log) though, this should be carefully benchmarked
 * though as my guess is that the improvement is going to be marginal and will depend on the
 * architecture (how well pow, exp and log can be optimized on a given architecture).
 *
 * If you have any experience with pow please suggest any better solution and/or provide a
 * faster implementation.
 */

struct InversePower_interaction {
    double const _eps;
    double const _pow;
    Array<double> const _radii;

    InversePower_interaction(double pow, double eps, Array<double> const radii)
    : _eps(eps),
      _pow(pow),
      _radii(radii.copy())
    {}

    /* calculate energy from distance squared */
    double energy(double r2, size_t atomi, size_t atomj) const
    {
        double E;
        double r0 = _radii[atomi] + _radii[atomj]; //sum of the hard core radii
        double r = std::sqrt(r2);
        if (r >= r0){
            E = 0.;
        }
        else{
            E = std::pow((1 -r/r0), _pow) * _eps/_pow;
        }
        return E;
    }

    /* calculate energy and gradient from distance squared, gradient is in -(dv/drij)/|rij| */
    double energy_gradient(double r2, double *gij, size_t atomi, size_t atomj) const
    {
        double E;
        double r0 = _radii[atomi] + _radii[atomj]; //sum of the hard core radii
        double r = std::sqrt(r2);

        if (r >= r0){
            E = 0.;
            *gij = 0;
        }
        else {
            double factor = std::pow((1 -r/r0), _pow) * _eps;
            E =  factor / _pow;
            *gij =  - factor / ((r-r0)*r);
        }

        return E;
    }

    double inline energy_gradient_hessian(double r2, double *gij, double *hij, size_t atomi, size_t atomj) const
    {
        double E;
        double r0 = _radii[atomi] + _radii[atomj]; //sum of the hard core radii
        double r = std::sqrt(r2);
        if (r >= r0){
            E = 0.;
            *gij = 0;
            *hij=0;
        }
        else {
            double factor = std::pow((1 -r/r0), _pow) * _eps;
            double denom = 1.0 / (r-r0);
            E =  factor / _pow;
            *gij =  - factor * denom / r ;
            *hij = (_pow-1) * factor * denom * denom;
        }

        return E;
    }
};


//
// combine the components (interaction, looping method, distance function) into
// defined classes
//

/**
 * Pairwise Inverse Power potential
 */
template <size_t ndim>
class InversePower : public SimplePairwisePotential< InversePower_interaction, cartesian_distance<ndim>> {
public:
    InversePower(double pow, double eps, pele::Array<double> const radii)
        : SimplePairwisePotential< InversePower_interaction, cartesian_distance<ndim>>(
                std::make_shared<InversePower_interaction>(pow, eps, radii),
                std::make_shared<cartesian_distance<ndim>>()
          )
          {}
};

template <size_t ndim>
class InversePowerPeriodic : public SimplePairwisePotential< InversePower_interaction, periodic_distance<ndim>> {
public:
    InversePowerPeriodic(double pow, double eps, pele::Array<double> const radii, pele::Array<double> const boxvec)
        : SimplePairwisePotential< InversePower_interaction, periodic_distance<ndim>>(
                std::make_shared<InversePower_interaction>(pow, eps, radii),
                std::make_shared<periodic_distance<ndim>>(boxvec)
          )
          {}
};

}
#endif


