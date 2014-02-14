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

namespace pele{

/*Random coords displacement, generates a random displacement for a N dimensional system
 * sampling from a N-dimensional sphere
 * */

class RandomCoordsDisplacement:public TakeStep{
protected:
	size_t _seed;
	std::mt19937_64 _generator;
	std::normal_distribution<double> _distribution;
	pele::Array<double> _dr;
	size_t _N;
public:
	RandomCoordsDisplacement(size_t N);
	virtual ~RandomCoordsDisplacement() {}
	virtual void takestep(Array<double>& coords, double stepsize, MC * mc);
	virtual size_t get_seed(){return _seed;}
};

RandomCoordsDisplacement::RandomCoordsDisplacement(size_t N):
		_seed(std::chrono::system_clock::now().time_since_epoch().count()),
		_generator(_seed), _distribution(0.0,1.0), _dr(N),
		_N(N)
		{
			std::cout<<"seed TakeStep:"<<_seed<<std::endl;
		}

void RandomCoordsDisplacement::takestep(Array<double>& coords, double stepsize, MC * mc){
	double norm2, rand;
	assert(coords.size() == _N);
	for(size_t i=0; i<_N;++i){
		rand = _distribution(_generator);
		_dr[i] = rand;
		norm2 += rand*rand;
	}

	_dr *= stepsize/sqrt(norm2);
	coords += _dr;
}

}
#endif
