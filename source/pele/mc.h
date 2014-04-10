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
 * _nitercount is the cumulative number of MCMC steps taken by the class
 * _neval is the number of energy evaluations
 * _stepsize the is the stepsize to pass to takestep
 * _temperature is the temperature at which the simulation is performed
 * _energy is the current energy of the system
 * _success records whether the step has been accepted or rejected
 */

class MC {
protected:
	Array<double> _coords, _trial_coords;
	shared_ptr<pele::BasePotential> _potential;
	list< shared_ptr<Action> > _actions;
	list< shared_ptr<AcceptTest> > _accept_tests;
	list< shared_ptr<ConfTest> > _conf_tests;
	shared_ptr<TakeStep> _takestep;
	size_t _niter, _nitercount, _accept_count, _E_reject_count, _conf_reject_count;
	bool _success;
	/*nitercount is the cumulative count, it does not get reset at the end of run*/
public:
	/*need to keep these public to make them accessible to tests and actions*/
	size_t _neval;
	double _stepsize, _temperature, _energy, _trial_energy;

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
	double get_trial_energy(){return _trial_energy;}
	Array<double> get_coords(){return _coords;}
	double get_accepted_fraction(){return ((double) _accept_count)/(_accept_count+_E_reject_count+_conf_reject_count);};
	double get_conf_rejection_fraction(){return ((double)_conf_reject_count)/_nitercount;};
	size_t get_iterations_count(){return _nitercount;};
	size_t get_neval(){return _neval;};
	double get_stepsize(){return _stepsize;};
};

MC::MC(pele::BasePotential * potential, Array<double>& coords, double temperature, double stepsize):
		_coords(coords.copy()),_trial_coords(_coords.copy()), _potential(potential),
			_niter(0), _nitercount(0), _accept_count(0), _E_reject_count(0),
			_conf_reject_count(0), _success(true), _neval(0), _stepsize(stepsize), _temperature(temperature)

		{
			_energy = _potential->get_energy(_coords);
			_trial_energy = _energy;
			//std::cout<<"Energy is "<<_energy<<std::endl;
			++_neval;
		}

void MC::one_iteration()
{
	_success = true;

	for(size_t i=0; i<_coords.size();++i){
		_trial_coords[i] = _coords[i];
	}

	_takestep->takestep(_trial_coords, _stepsize, this);

	for (auto& test : _conf_tests ){
		_success = test->test(_trial_coords, this);
		if (_success == false)
			++_conf_reject_count;
			break;
	}

	if (_success == true)
	{
		_trial_energy = _potential->get_energy(_trial_coords);
		++_neval;
		
		for (auto& test : _accept_tests ){
			_success = test->test(_trial_coords, _trial_energy, _coords, _energy, _temperature, this);
			if (_success == false)
			{
				++_E_reject_count;
				break;
			}
		}

		if (_success == true){
			for(size_t i=0;i<_coords.size();++i)
			{
			  _coords[i] = _trial_coords[i];
			}
			_energy = _trial_energy;
			++_accept_count;
		}

		/*currently if conf test is rejected we don't record the state*/
		for (auto& action : _actions){
				action->action(_coords, _energy, _success, this);
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
