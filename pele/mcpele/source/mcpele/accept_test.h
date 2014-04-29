#ifndef _PELE_ACCEPT_TEST_H__
#define _PELE_ACCEPT_TEST_H__

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
	MetropolisTest(size_t rseed);
	virtual ~MetropolisTest() {}
	virtual bool test(Array<double> &trial_coords, double trial_energy, Array<double> & old_coords, double old_energy, double temperature, MC * mc);
	virtual size_t get_seed(){return _seed;}
};

MetropolisTest::MetropolisTest(size_t rseed):
		_seed(rseed), _generator(_seed), _distribution(0.0,1.0)
		{
        #ifdef DEBUG
			std::cout<<"seed Metropolis:"<<_seed<<std::endl;
			//std::chrono::system_clock::now().time_since_epoch().count()
        #endif
		}

bool MetropolisTest::test(Array<double> &trial_coords, double trial_energy, Array<double>& old_coords, double old_energy, double temperature, MC * mc)
{
	double rand, w, wcomp;
	bool success = true;

	wcomp = (trial_energy - old_energy) / temperature;
	w = exp(-wcomp);

	if (w < 1.0)
	{
		rand = _distribution(_generator);
		if (rand > w)
			success = false;
	}

	return success;
}

/*ENERGY WINDOW TEST
 * this test checks that the energy of the system stays within a certain energy window
 * */

class EnergyWindowTest:public AcceptTest{
protected:
	double _min_energy, _max_energy;
public:
	EnergyWindowTest(double min_energy, double max_energy);
	virtual ~EnergyWindowTest() {}
	virtual bool test(Array<double> &trial_coords, double trial_energy, Array<double> & old_coords, double old_energy, double temperature, MC * mc);
};

EnergyWindowTest::EnergyWindowTest(double min_energy, double max_energy):
		_min_energy(min_energy),_max_energy(max_energy)
		{}

bool EnergyWindowTest::test(Array<double> &trial_coords, double trial_energy, Array<double> & old_coords, double old_energy, double temperature, MC * mc)
{
	bool success = true;

	if ((trial_energy < _min_energy) or (trial_energy > _max_energy))
		success = false;

	return success;
}

}
#endif
