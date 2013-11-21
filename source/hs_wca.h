#ifndef _PELE_HS_WCA_H
#define _PELE_HS_WCA_H

#include "pairwise_potential.h"
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
        double energy(double r2, size_t i, size_t j) const {
            double E;
        	double r = sqrt(r2);
        	double r0 = _radii[i] + _radii[j]; //sum of the hard core radii
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
            	std::cout<<"WARNING: distance between particles "<<i<<" and "<<j<<" is less than their hard core radii"<<std::endl;
            }
            else if(r < coff )
            	E = _eps*(-C6*ir6 + C12*ir12 + 1.0/4);
            else
            	E = 0.;

            return E;
        }

        /* calculate energy and gradient from distance squared, gradient is in g/|rij|, r0 is the hard core distance, r is the distance between the centres */
        double energy_gradient(double r2, double *gij, size_t i, size_t j) const {
        	double E;
			double r = sqrt(r2);
			double r0 = _radii[i] + _radii[j]; //sum of the hard core radii
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
				std::cout<<"WARNING: distance between particles "<<i<<" and "<<j<<" is less than their hard core radii"<<std::endl;
			}
			else if(r < coff)
			{
				E = _eps*(- C6 * ir6 + C12 * ir12 + 1.0/4);
				*gij = _eps*(12 * C12 * ir12 - 6 * C6 * ir6) * ir2;
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
    class HS_WCA : public PairwisePotential< HS_WCA_interaction >
    {
        public:
            HS_WCA(double eps, double sca, Array<double> radii)
                : PairwisePotential< HS_WCA_interaction > ( new HS_WCA_interaction(eps, sca, radii) ) {}
    };

    /**
     * Pairwise HS_WCA potential in a rectangular box
     */
    class HS_WCAPeriodic : public PairwisePotential< HS_WCA_interaction, periodic_distance > {
        public:
            HS_WCAPeriodic(double eps, double sca, Array<double> radii, double const *boxvec)
                : PairwisePotential< HS_WCA_interaction, periodic_distance> (
                        new HS_WCA_interaction(eps, sca, radii),
                        new periodic_distance(boxvec[0], boxvec[1], boxvec[2])
                        )
            {}
    };
}
#endif
