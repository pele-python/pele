#ifndef _PELE_MC_H
#define _PELE_MC_H

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
using std::sqrt;

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
 * Accept Test
 */

class ConfTest{
public:
	virtual ~ConfTest(){}
	virtual bool test(Array<double> &trial_coords, MC * mc) =0;
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
 * _coords and _trialcoords are arrays that store coordinates and trial coordinates respectively
 * _potential is on object of Pele::BasePotential type that defines the interaction potential
 * _E_reject_count is the count of rejections due to an energy test (e.g. Metropolis)
 * _conf_reject_count is the count of rejections due to a configuration test (e.g. spherical container)
 * _niter is the count of steps whithin a MCMC run, it is reset to zero at the end of the run
 * _nitercoutn is the cumulative number of MCMC steps taken by the class
 * _neval is the number of energy evaluations
 * _stepsize the is the stepsize to pass to takestep
 * _temperature is the temperature at which the simulation is performed
 * _energy is the current energy of the system
 */

class MC {
protected:
	Array<double> _coords, _trial_coords;
	shared_ptr<pele::BasePotential> _potential;
	list< shared_ptr<Action> > _actions;
	list< shared_ptr<AcceptTest> > _accept_tests;
	list< shared_ptr<ConfTest> > _conf_tests;
	shared_ptr<TakeStep> _takestep;
	size_t _niter, _nitercount, _neval, _accept_count, _E_reject_count, _conf_reject_count;
	/*nitercount is the cumulative count, it does not get reset at the end of run*/
public:
	/*need to keep these public to make them accessible to tests and actions*/
	double _stepsize, _temperature, _energy;

	MC(pele::BasePotential * potential, Array<double>& coords, double temperature, double stepsize);

	~MC(){}

	void one_iteration();
	void run(size_t max_iter);
	void set_temperature(double T){_temperature = T;}
	void set_stepsize(double stepsize){_stepsize = stepsize;}
	void add_action(shared_ptr<Action> action){_actions.push_back(action);}
	void add_accept_test( shared_ptr<AcceptTest> accept_test){_accept_tests.push_back(accept_test);}
	void add_conf_test( shared_ptr<ConfTest> conf_test){_conf_tests.push_back(conf_test);}
	void set_takestep( shared_ptr<TakeStep> takestep){_takestep = takestep;}
	void set_coordinates(Array<double>& coords, double energy){
		_coords = coords.copy();
		_energy = energy;
	}
	double get_energy(){return _energy;}
	Array<double> get_coords(){return _coords;}
	double get_accepted_fraction(){return ((double) _accept_count)/(_accept_count+_E_reject_count+_conf_reject_count);};
	size_t get_iterations_count(){return _nitercount;};
};

MC::MC(pele::BasePotential * potential, Array<double>& coords, double temperature, double stepsize):
		_coords(coords.copy()),_trial_coords(_coords.copy()), _potential(potential),
			_niter(0), _nitercount(0), _neval(0), _accept_count(0),_E_reject_count(0),
			_conf_reject_count(0), _stepsize(stepsize), _temperature(temperature)

		{
			_energy = _potential->get_energy(_coords);
			//std::cout<<"Energy is "<<_energy<<std::endl;
			++_neval;
		}

void MC::one_iteration()
{
	double trial_energy;
	bool success = true;

	for(size_t i=0; i<_coords.size();++i){
		_trial_coords[i] = _coords[i];
	}

	_takestep->takestep(_trial_coords, _stepsize, this);

	for (auto& test : _conf_tests ){
		success = test->test(_trial_coords, this);
		if (success == false)
			++_conf_reject_count;
			break;
	}

	if (success == true)
	{
		trial_energy = _potential->get_energy(_trial_coords);
		++_neval;
		
		for (auto& test : _accept_tests ){
			success = test->test(_trial_coords, trial_energy, _coords, _energy, _temperature, this);
			if (success == false)
			{
				++_E_reject_count;
				break;
			}
		}

		if (success == true){
			for(size_t i=0;i<_coords.size();++i)
			{
			  _coords[i] = _trial_coords[i];
			}
			_energy = trial_energy;
			++_accept_count;
		}

		/*consider moving actions outside the accepted configuration test condition?
		(currently if conf test is rejected we don't record the state)*/
		for (auto& action : _actions){
				action->action(_coords, _energy, success, this);
			}
	}

	++_niter;
	++_nitercount;
}

void MC::run(size_t max_iter)
{
	while(_niter < max_iter)
		this->one_iteration();
	_niter = 0;
}

}

#endif
