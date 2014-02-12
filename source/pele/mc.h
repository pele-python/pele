#ifndef _PELE_MC_H__
#define _PELE_MC_H__

#include <math.h>
#include <algorithm>
#include <list>
#include "array.h"
#include "base_potential.h"

using std::vector;
using std::list;
using std::runtime_error;
using pele::Array;

namespace pele{

class MC;

/*
 * Action
 */

class Action {
public:
	virtual ~Action(){}
	virtual void action(Array<double> &coords, double energy, bool accepted, MC * mc) =0;
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
	Array<double> _coords;
	pele::BasePotential * _potential;
	list<pele::Action*> _actions;
	list<pele::AcceptTest*> _accept_tests;
	pele::TakeStep * _takestep;
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
	void add_action(Action* action){_actions.push_back(action);}
	void add_accept_test(AcceptTest* accept_test){_accept_tests.push_back(accept_test);}
	void set_coordinates(Array<double> coords, double energy){
		_coords = coords;
		_energy = energy;
	}
	Array<double> get_coordinates(){return _coords;}
	double get_energy(){return _energy;}

};

MC::MC(pele::BasePotential * potential, Array<double> coords, double temperature, double stepsize):
		_stepsize(stepsize),_temperature(temperature), _niter(0), _neval(0),
		_coords(coords),_potential(potential)
		{
			_energy = _potential->get_energy(_coords);
			++_neval;
		}
void MC::one_iteration()
{
	Array<double> trial_coords;
	double trial_energy;
	bool success = true;

	trial_coords = _coords.copy();

	_takestep->takestep(trial_coords, _stepsize, this);

	trial_energy = _potential->get_energy(trial_coords);
	++_neval;

	for (list<AcceptTest*>::iterator test = _accept_tests.begin(); test != _accept_tests.end(); ++test){

		success *= (*test)->test(trial_coords, trial_energy, _coords, _energy, _temperature, this);
		if (success == false)
			break;
	}

	if (success == true){
		_coords = trial_coords.copy();
		_energy = trial_energy;
	}

	for (list<Action*>::iterator action = _actions.begin(); action != _actions.end(); ++action){

			(*action)->action(_coords, _energy, success, this);
	}

	++_niter;
}

void MC::run(size_t max_iter)
{
	while(_niter < max_iter)
		this->one_iteration();
}

/*
 *
 * Inherited classes (specific actions, tests, takesteps)
 *
 * /


/*Adjust Step*/

class AdjustStep : public Action {
public:
	double _factor;
	AdjustStep(double factor);
	virtual ~AdjustStep() {}
	virtual void action(Array<double> coords, double energy, bool accepted, MC * mc);
};

AdjustStep::AdjustStep(double factor):
			_factor(factor){}

void AdjustStep::action(Array<double> coords, double energy, bool accepted, MC * mc) {
		if (accepted == false)
			mc->_stepsize *= _factor;
		else
			mc->_stepsize /= _factor;
	}

}
#endif
