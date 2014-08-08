#ifndef _PELE_INVERSEPOWER_H
#define _PELE_INVERSEPOWER_H

#include "simple_pairwise_potential.h"
#include "simple_pairwise_ilist.h"
#include "atomlist_potential.h"
#include "distance.h"
#include <memory>

namespace pele {

/**
 * Pairwise interaction for Weeks-Chandler-Andersen (WCA) potential
 */

struct InversePower_interaction {
    double const _eps; //cutoff distance for WCA potential
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
        if (r >= r0)
            E = 0.;
        else
        	E = std::pow((1 -r/r0), _pow) * _eps/_pow;

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
			*gij =  - factor / (r-r0);
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
			*gij =  - factor * denom ;
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


