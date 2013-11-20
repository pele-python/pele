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
     */
    struct HS_WCA_interaction {
    	double const _eps, _sca;
    	double const _infty, _prfac;

        HS_WCA_interaction(double eps, double sca) :
            _eps(eps), _sca(1 + 1./sca),
            _infty(pow(10.0,10)), _prfac(1./(sqrt(2)*sca*sca*sca))
    	{}
        // _prfac is the cubic power of 1/(2**(1/6) * sca)
        //well depth _eps and scaling factor (shell thickness = 1/sca * R, where R is the hard core radius)

        /* calculate energy from distance squared, r0 is the hard core distance, r is the distance between the centres */
        double energy(double r2, double r0) const {
            double E;
        	double r = sqrt(r2);
        	double dr = r - r0;
            double ir2 = 1.0/(dr*dr);
            double ir6 = ir2*ir2*ir2;
            double ir12 = ir6*ir6;
            double C3 = _prfac*r0*r0*r0;
            double C6 = C3*C3;
            double C12 = C6*C6;
            double coff = r0*_sca; //distance at which the soft cores are at contact
            if (r <= r0)
            	E = _infty;
            else if(r < coff )
            	E = _eps*(-C6*ir6 + C12*ir12 + 1.0/4);
            else
            	E = 0.;

            return E;
        }

        /* calculate energy and gradient from distance squared, gradient is in g/|rij|, r0 is the hard core distance, r is the distance between the centres */
        double energy_gradient(double r2, double r0, double *gij) const {
        	double E;
			double r = sqrt(r2);
			double dr = r - r0;
			double ir2 = 1.0/(dr*dr);
			double ir6 = ir2*ir2*ir2;
			double ir12 = ir6*ir6;
			double C3 = _prfac*r0*r0*r0;
			double C6 = C3*C3;
			double C12 = C6*C6;
			double coff = r0*_sca; //distance at which the soft cores are at contact

			if (r <= r0)
			{
				E = _infty;
				*gij = _infty;
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
    class HS_WCA : public SimplePairwisePotential< HS_WCA_interaction >
    {
        public:
            HS_WCA(double eps, double sca)
                : SimplePairwisePotential< HS_WCA_interaction > ( new HS_WCA_interaction(eps, sca) ) {}
    };

    /**
     * Pairwise HS_WCA potential in a rectangular box
     */
    class HS_WCAPeriodic : public SimplePairwisePotential< HS_WCA_interaction, periodic_distance > {
        public:
            HS_WCAPeriodic(double eps, double sca, double const *boxvec)
                : SimplePairwisePotential< HS_WCA_interaction, periodic_distance> (
                        new HS_WCA_interaction(eps, sca),
                        new periodic_distance(boxvec[0], boxvec[1], boxvec[2])
                        )
            {}
    };

    /**
     * Pairwise HS_WCA potential with interaction lists
     */
    class HS_WCA_interaction_list : public SimplePairwiseInteractionList< HS_WCA_interaction > {
        public:
            HS_WCA_interaction_list(Array<long int> & ilist, double eps, double sca)
                :  SimplePairwiseInteractionList< HS_WCA_interaction > ( new HS_WCA_interaction(eps, sca), ilist) {}
    };
}
#endif
