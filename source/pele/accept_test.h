#ifndef _ACCEPT_TEST_H__
#define _ACCEPT_TEST_H__

#include <math.h>
#include <algorithm>
#include <random>
#include <chrono>
#include "array.h"
#include "mc.h"


using std::runtime_error;
using pele::Array;

namespace pele{

/*Metropolis acceptance criterion
 * */

class MetropolisTest:public AcceptTest{
protected:
	size_t _seed;
	std::mt19937_64 _generator;
	std::uniform_real_distribution<double> _distribution;
public:
	MetropolisTest();
	virtual ~MetropolisTest() {}
	virtual bool test(Array<double> &trial_coords, double trial_energy, Array<double> & old_coords, double old_energy, double temperature, MC * mc);
	virtual size_t get_seed(){return _seed;}
};

MetropolisTest::MetropolisTest():
		_seed(std::chrono::system_clock::now().time_since_epoch().count()),
		_generator(_seed), _distribution(0.0,1.0)
		{}

bool MetropolisTest::test(Array<double> &trial_coords, double trial_energy, Array<double> & old_coords, double old_energy, double temperature, MC * mc)
{
	double rand, w, wcomp;
	bool success = true;

	wcomp = (trial_energy - old_energy) / temperature;
	w = std::min(1.0,exp(-wcomp));

	if (w < 1.0)
	{
		rand = _distribution(_generator);
		if (rand > w)
			success = false;
	}

	return success;
}

}
#endif
