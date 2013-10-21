#ifndef PYGMIN_SIMPLE_PAIRWISE_POTENTIAL_H
#define PYGMIN_SIMPLE_PAIRWISE_POTENTIAL_H

#include "potential.h"

namespace pele
{
	template<typename pairwise_interaction>
	class SimplePairwisePotential : public Potential
	{
	protected:
		pairwise_interaction *_interaction;
		SimplePairwisePotential(pairwise_interaction *interaction) : _interaction(interaction) {}

	public:
		virtual double get_energy(Array x);
		virtual double get_energy_gradient(Array x, Array grad);
	};

	template<typename pairwise_interaction>
	inline double SimplePairwisePotential<pairwise_interaction>::get_energy_gradient(Array x, Array grad)
	{
		double e=0.;
		double gij, dr[3];
		size_t const natoms = x.size()/3;

		for(size_t i=0; i<x.size(); ++i)
			grad[i] = 0.;

		for(int i=0; i<natoms; ++i) {
			int i1 = 3*i;
			for(int j=i+1; j<natoms; ++j) {
				int i2 = 3*j;
				for(int k=0; k<3; ++k)
					dr[k] = x(i1+k) - x(i2+k);

				double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];

				e += _interaction->energy_gradient(r2, &gij);
				for(size_t k=0; k<3; ++k)
					grad[i1+k] -= gij * dr[k];
				for(size_t k=0; k<3; ++k)
					grad[i2+k] += gij * dr[k];
			}
		}

		return e;
	}

	template<typename pairwise_interaction>
	inline double SimplePairwisePotential<pairwise_interaction>::get_energy(Array x)
	{
		double e=0.;
		size_t const natoms = x.size()/3;

		for(size_t i=0; i<natoms; ++i) {
			size_t i1 = 3*i;
			for(size_t j=i+1; j<natoms; ++j) {
				size_t i2 = 3*j;
				double dr[3];
				for(size_t k=0; k<3; ++k)
					dr[k] = x(i1+k) - x(i2+k);
				double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
				e += _interaction->energy(r2);
			}
		}

		return e;
	}
}

#endif
