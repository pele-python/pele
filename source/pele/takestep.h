#ifndef _PELE_TAKESTEP_H__
#define _PELE_TAKESTEP_H__

#include <math.h>
#include <algorithm>
#include <random>
#include <chrono>
#include "array.h"
#include "mc.h"


using std::runtime_error;
using pele::Array;
using std::sqrt;

namespace pele{

/*Random coords displacement, generates a random displacement for a N dimensional system
 * sampling from a N-dimensional sphere
 * the stepsize is defined per coordinates, that's why the maximum stepsize is sqrt(N)*stepsize
 * */

class RandomCoordsDisplacement:public TakeStep{
protected:
	size_t _seed;
	std::mt19937_64 _generator;
	std::uniform_real_distribution<double> _distribution;
	size_t _N;
	double _rootN;
public:
	RandomCoordsDisplacement(size_t N, size_t rseed);
	virtual ~RandomCoordsDisplacement() {}
	virtual void takestep(Array<double>& coords, double stepsize, MC * mc);
	virtual size_t get_seed(){return _seed;}
};

RandomCoordsDisplacement::RandomCoordsDisplacement(size_t N, size_t rseed):
		_seed(rseed),_generator(_seed), _distribution(0.0,1.0), _N(N)
		{
        #ifdef DEBUG
			std::cout<<"seed TakeStep:"<<_seed<<std::endl;
        #endif
		}

void RandomCoordsDisplacement::takestep(Array<double>& coords, double stepsize, MC * mc){
	double norm2=0;
	double rand;
	//assert(coords.size() == _N);
	for(size_t i=0; i<_N;++i){
		rand = _distribution(_generator);
		coords[i] += (0.5-rand)*stepsize;
	}
}

}
#endif
