#ifndef _PELE_HARMONIC_H
#define _PELE_HARMONIC_H

#include "base_potential.h"
#include "distance.h"

namespace pele {

template<typename distance_policy = cartesian_distance >
    class BaseHarmonic : public BasePotential {
    protected:
    	distance_policy *_dist;
    	pele::Array<double> _origin;
    	double _k;
    	BaseHarmonic(pele::Array<double> origin, double k, distance_policy *dist=NULL):
		_dist(dist), _origin(origin.copy()), _k(k)
			{
				if(_dist == NULL) _dist = new distance_policy;
			}
    public:
    	~BaseHarmonic(){
    		if (_dist != NULL) delete _dist;
    	}
    	double get_energy(pele::Array<double> x);
    	double get_energy_gradient(pele::Array<double> x, pele::Array<double> grad);
    	void set_k(double newk){_k = newk;};
    };


    /* calculate energy from distance squared, r0 is the hard core distance, r is the distance between the centres */
    template<typename distance_policy>
    inline double BaseHarmonic<distance_policy>::get_energy(pele::Array<double> x) {
    	assert(x.size() == _origin.size());
    	double norm2 = 0;
		double dr[3];
		size_t N = x.size()/3;
		for(size_t i=0; i<N;++i)
		{
			size_t i1 = i*3;
			_dist->get_rij(dr, &x[i1], &_origin[i1]);
			norm2 += dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
		}

		return 0.5 * _k * norm2;
	}

	/* calculate energy and gradient from distance squared, gradient is in g/|rij|, r0 is the hard core distance, r is the distance between the centres */
    template<typename distance_policy>
    inline double BaseHarmonic<distance_policy>::get_energy_gradient(pele::Array<double> x, pele::Array<double> grad) {
    	assert(x.size() == _origin.size());
    	assert(grad.size() == _origin.size());
    	double norm2 = 0;
		double dr[3];
		size_t N = x.size()/3;

		for(size_t i=0; i<N;++i)
		{
			size_t i1 = i*3;
			_dist->get_rij(dr, &x[i1], &_origin[i1]);
			norm2 += dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
			for(size_t k=0; k<3; ++k)
				grad[i1+k] = _k * dr[k];
		}

		return 0.5 * _k * norm2;
	}

    /**
	 * Pairwise HS_WCA potential
	 */
	class Harmonic : public BaseHarmonic<>
	{
		public:
			Harmonic(pele::Array<double> origin, double k)
				: BaseHarmonic(origin, k){}
	};

	/**
	 * Pairwise HS_WCA potential in a rectangular box
	 */
	class HarmonicPeriodic : public BaseHarmonic<periodic_distance> {
		public:
		HarmonicPeriodic(pele::Array<double> origin, double k, double const *boxvec)
				: BaseHarmonic<periodic_distance> (
						origin, k,
						new periodic_distance(boxvec[0], boxvec[1], boxvec[2])
						)
			{}
	};


}

#endif
