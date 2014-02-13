#ifndef _PELE_MC_H__
#define _PELE_MC_H__

#include <math.h>
#include <algorithm>
#include <list>
#include <memory>
#include "array.h"
#include "base_potential.h"

using std::list;
using std::runtime_error;
using std::shared_ptr;
using pele::Array;

namespace pele{

class MC;

/*
 * Action
 */

class Action {
public:
	virtual ~Action(){}
	virtual void action(Array<double> &coords, double energy, bool accepted, MC* mc) =0;
};

/*
 * Accept Test
 */

class AcceptTest{
public:
	virtual ~AcceptTest(){}
	virtual bool test(Array<double> &trial_coords, double trial_energy, Array<double> & old_coords, double old_energy, double temperature, MC * mc) =0;
};

/*
 * Take Step
 */

class TakeStep{
public:
	virtual ~TakeStep(){}
	virtual void takestep(Array<double> &coords, double stepsize, MC * mc) =0;
};

/*
 * Monte Carlo
 */

class MC {
protected:
	Array<double> _coords, _trial_coords;
	shared_ptr<pele::BasePotential> _potential;
	list< shared_ptr<Action> > _actions;
	list< shared_ptr<AcceptTest> > _accept_tests;
	list< shared_ptr<AcceptTest> > _conf_tests;
	shared_ptr<TakeStep> _takestep;
	size_t _niter, _neval;
public:
	/*need to keep these public to make them accessible to tests and actions*/
	double _stepsize, _temperature, _energy;

	MC(pele::BasePotential * potential, Array<double> coords, double temperature, double stepsize);

	~MC(){}

	void one_iteration();
	void run(size_t max_iter);
	void set_temperature(double T){_temperature = T;}
	void set_stepsize(double stepsize){_stepsize = stepsize;}
	void add_action(shared_ptr<Action> action){_actions.push_back(action);}
	void add_accept_test( shared_ptr<AcceptTest> accept_test){_accept_tests.push_back(accept_test);}
	void add_conf_test( shared_ptr<AcceptTest> conf_test){_conf_tests.push_back(conf_test);}
	void set_take_step( shared_ptr<TakeStep> takestep){_takestep = takestep;}
	void set_coordinates(Array<double> coords, double energy){
		_coords = coords;
		_energy = energy;
	}
	Array<double> get_coordinates(){return _coords;}
	double get_energy(){return _energy;}

};

MC::MC(pele::BasePotential * potential, Array<double> coords, double temperature, double stepsize):
		_stepsize(stepsize),_temperature(temperature), _niter(0), _neval(0),
		_coords(coords), _trial_coords(_coords.size()), _potential(potential)
		{
			_energy = _potential->get_energy(_coords);
			++_neval;
		}

void MC::one_iteration()
{
	Array<double> trial_coords;
	double trial_energy;
	bool success = true;

	for(int i=0;_coords.size();++i){
		_trial_coords[i] = _coords[i];
	}

	_takestep->takestep(_trial_coords, _stepsize, this);

	for (auto& test : _conf_tests ){

			success = test->test(_trial_coords, trial_energy, _coords, _energy, _temperature, this);
			if (success == false)
				break;
		}

	trial_energy = _potential->get_energy(_trial_coords);
	++_neval;

	for (auto& test : _accept_tests ){

		success = test->test(_trial_coords, trial_energy, _coords, _energy, _temperature, this);
		if (success == false)
			break;
	}

	if (success == true){
		for(int i=0;_coords.size();++i)
		{
			_coords[i] = _trial_coords[i];
		}
		_energy = trial_energy;
	}

	for (auto& action : _actions){

			action->action(_coords, _energy, success, this);
	}

	++_niter;
}

void MC::run(size_t max_iter)
{
	while(_niter < max_iter)
		this->one_iteration();
}

}

#endif
