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
     * Pairwise interaction for Hard Sphere + Weeks-Chandler-Andersen (HS_WCA) potential, refer to D. Asenjo PhD thesis pp 66
     * _prfac is the cubic power of _sca/(2**(1/6))
     * well depth _eps and scaling factor (shell thickness = sca * R, where R is the hard core radius), sca determined the thickness of the shell
     */

    struct HS_WCA_interaction {
    	double const _eps, _sca;
    	double const _infty, _prfac;
    	Array<double> const _radii;

        HS_WCA_interaction(double eps, double sca, Array<double> radii) :
            _eps(eps), _sca(sca),
            _infty(pow(10.0,10)), _prfac(sca*sca*sca/sqrt(2)),
            _radii(radii.copy())
    	{}

        /* calculate energy from distance squared, r0 is the hard core distance, r is the distance between the centres */
        double energy(double r2, size_t atomi, size_t atomj) const {
            double E;
        	double r = sqrt(r2);
        	double r0 = _radii[atomi] + _radii[atomj]; //sum of the hard core radii
        	double dr = r - r0;
            double ir2 = 1.0/(dr*dr);
            double ir6 = ir2*ir2*ir2;
            double ir12 = ir6*ir6;
            double C3 = _prfac*r0*r0*r0;
            double C6 = C3*C3;
            double C12 = C6*C6;
            double coff = r0*(1.0 +_sca); //distance at which the soft cores are at contact
            if (r <= r0)
            {
            	E = _infty;
            	std::cout<<"WARNING: distance between atoms "<<atomi<<" and "<<atomj<<" is "<<r0-r<<", less than their hard core separation"<<std::endl;
            }
            else if(r < coff )
            	E = 4.*_eps*(-C6*ir6 + C12*ir12) + _eps;
            else
            	E = 0.;

            return E;
        }

        /* calculate energy and gradient from distance squared, gradient is in g/|rij|, r0 is the hard core distance, r is the distance between the centres */
        double energy_gradient(double r2, double *gij, size_t atomi, size_t atomj) const {
        	double E;
			double r = sqrt(r2);
			double r0 = _radii[atomi] + _radii[atomj]; //sum of the hard core radii
			double dr = r - r0;
			double ir2 = 1.0/(dr*dr);
			double ir6 = ir2*ir2*ir2;
			double ir12 = ir6*ir6;
			double C3 = _prfac*r0*r0*r0;
			double C6 = C3*C3;
			double C12 = C6*C6;
			double coff = r0*(1+_sca); //distance at which the soft cores are at contact

			if (r <= r0)
			{
				E = _infty;
				*gij = _infty;
				std::cout<<"WARNING: distance between atoms "<<atomi<<" and "<<atomj<<" is "<<r0-r<<"less than their hard core separation"<<std::endl;
			}
			else if(r < coff)
			{
				E = 4.*_eps*(- C6 * ir6 + C12 * ir12) + _eps;
				*gij = 4.*_eps*(- 6 * C6 * ir6 + 12 * C12 * ir12) / (dr*r); // this is -g|gij| (for consistency with the loop in pairwise potential)
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
            HS_WCA(double eps, double sca, Array<double> radii)
                : SimplePairwisePotential< HS_WCA_interaction > ( new HS_WCA_interaction(eps, sca, radii) ) {}
    };

    /**
     * Pairwise HS_WCA potential in a rectangular box
     */
    class HS_WCAPeriodic : public SimplePairwisePotential< HS_WCA_interaction, periodic_distance > {
        public:
            HS_WCAPeriodic(double eps, double sca, Array<double> radii, double const *boxvec)
                : SimplePairwisePotential< HS_WCA_interaction, periodic_distance> (
                        new HS_WCA_interaction(eps, sca, radii),
                        new periodic_distance(boxvec[0], boxvec[1], boxvec[2])
                        )
            {}
    };

    /**
	 * Pairwise WCA potential with interaction lists
	 */
	class HS_WCANeighborList : public SimplePairwiseNeighborList< HS_WCA_interaction > {
		public:
			HS_WCANeighborList(Array<long int> & ilist, double eps, double sca, Array<double> radii)
				:  SimplePairwiseNeighborList< HS_WCA_interaction > ( new HS_WCA_interaction(eps, sca, radii), ilist) {}
	};
}
#endif
