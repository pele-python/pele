#ifndef _PELE_HARMONIC_H
#define _PELE_HARMONIC_H

#include "base_potential.h"
#include "distance.h"

namespace pele {

    class Harmonic : public BasePotential {
    protected:
    	double const _k;
    	pele::Array<double> _origin;
    public:
    	Harmonic(pele::Array<double> origin, double k);
    	~Harmonic(){}
    	double get_energy(pele::Array<double> x);
    	double get_energy_gradient(pele::Array<double> x, pele::Array<double> grad);
    };


    Harmonic::Harmonic(pele::Array<double> origin, double k) :
    		_k(k), _origin(origin.copy())
    		{}

    /* calculate energy from distance squared, r0 is the hard core distance, r is the distance between the centres */
	double Harmonic::get_energy(pele::Array<double> x) {
		double dx;
		double norm2 = 0;
		assert(x.size() == _origin.size());
		for(size_t i=0; i<x.size();++i)
		{
			dx = x[i] - _origin[i];
			norm2 += dx*dx;
		}

		return 0.5 * _k * norm2;
	}

	/* calculate energy and gradient from distance squared, gradient is in g/|rij|, r0 is the hard core distance, r is the distance between the centres */
	double Harmonic::get_energy_gradient(pele::Array<double> x, pele::Array<double> grad) {
		double dx;
		double norm2 = 0;
		assert(x.size() == _origin.size());
		assert(grad.size() == _origin.size());
		for(size_t i=0; i<x.size();++i)
		{
			dx = x[i] - _origin[i];
			grad[i] = _k * dx;
			norm2 += dx*dx;
		}

		return 0.5 * _k * norm2;
	}
}

#endif
