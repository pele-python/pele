#ifndef _PELE_WCA_H
#define _PELE_WCA_H

#include "simple_pairwise_potential.h"
#include "simple_pairwise_ilist.h"
#include "atomlist_potential.h"
#include "distance.h"

namespace pele {

/**
 * Pairwise interaction for Weeks-Chandler-Andersen (WCA) potential
 */
template <double alpha>
struct InversePower_interaction {
    double const _eps; //cutoff distance for WCA potential
    Array<double> const _radii;

    InversePower_interaction(double eps, Array<double> radii)
        : _eps(eps),
          _radii(radii.copy())
    {}

    /* calculate energy from distance squared */
    double energy(double r2, size_t atom_i, size_t atom_j) const 
    {
        double E;
        double r0 = _radii[atomi] + _radii[atomj]; //sum of the hard core radii
        double r = sqrt(r2);
        if (r < r0)
            E = std::pow((r/r0 - 1), alpha) * _eps/alpha;
        else
            E = 0.;

        return E;
    }

    /* calculate energy and gradient from distance squared, gradient is in -(dv/drij)/|rij| */
    double energy_gradient(double r2, double *gij, size_t atom_i, size_t atom_j) const 
    {
    	double E;
		double r0 = _radii[atomi] + _radii[atomj]; //sum of the hard core radii
		double r = sqrt(r2);
		if (r < r0){
			double factor = std::pow((r/r0 - 1), alpha) * _eps;
			E =  factor / alpha;
            *gij =  factor / (r-r0);
        }
		else {
            E = 0.;
            *gij = 0;
        }

        return E;
    }

    double inline energy_gradient_hessian(double r2, double *gij, double *hij, size_t atom_i, size_t atom_j) const 
    {
    	double E;
		double r0 = _radii[atomi] + _radii[atomj]; //sum of the hard core radii
		double r = sqrt(r2);
		if (r < r0){
			double factor = std::pow((r/r0 - 1), alpha) * _eps;
			double denom = 1.0 / (r-r0);
			E =  factor / alpha;
			*gij =  factor * denom ;
            *hij = (alpha-1) * factor * denom * denom;
        }
		else {
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


