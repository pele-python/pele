#ifndef _PELE_ACTIONS_H
#define _PELE_ACTIONS_H

#include <math.h>
#include <algorithm>
#include <list>
#include "array.h"
#include "mc.h"
#include "histogram.h"
#include "distance.h"

using std::runtime_error;
using pele::Array;
using std::sqrt;

namespace pele{

/*Adjust Step
 * 	factor is a multiplicative factor by which the stepsize is adjusted
 * 	niter determines the number of steps for which the action should take effect (generally
 * 	we want to adjust the step size only at the beginning of a simulation)
 * 	navg is the number of steps over which the acceptance is averaged
 * 	factor must be 0<f<1, if rejected make step shorter, if accepted make step longer
*/

class AdjustStep : public Action {
protected:
	double _target, _factor, _acceptedf;
	size_t _niter, _navg, _count, _naccepted, _nrejected;
public:
	AdjustStep(double target, double factor, size_t niter, size_t navg);
	virtual ~AdjustStep() {}
	virtual void action(Array<double> &coords, double energy, bool accepted, MC* mc);
};

AdjustStep::AdjustStep(double target, double factor, size_t niter, size_t navg):
			_target(target),_factor(factor),_acceptedf(0),
			_niter(niter),_navg(navg),_count(0),_naccepted(0),
			_nrejected(0){}


void AdjustStep::action(Array<double> &coords, double energy, bool accepted, MC* mc) {

	_count = mc->get_iterations_count();

	if (_count < _niter)
		{
			if (accepted == true)
				++_naccepted;
			else
				++_nrejected;

			if(_count % _navg == 0)
			{
				_acceptedf = (double) _naccepted / (_naccepted + _nrejected);

				//std::cout<<"acceptance "<<_acceptedf<<std::endl;
				//std::cout<<"stepsize before"<<mc->_stepsize<<std::endl;
				if (_acceptedf < _target)
					mc->_stepsize *= _factor;
				else
					mc->_stepsize /= _factor;
				//std::cout<<"stepsize after"<<mc->_stepsize<<std::endl;

				//now reset to zero memory of acceptance and rejection
				_naccepted = 0;
				_nrejected = 0;
			}

		}
}

/*Record energy histogram
*/

class RecordEnergyHistogram : public Action {
protected:
	pele::Histogram * _hist;
	double _bin;
	size_t _eqsteps, _count;
public:
	RecordEnergyHistogram(double min, double max, double bin, size_t eqsteps);
	virtual ~RecordEnergyHistogram() {delete _hist;}

	virtual void action(Array<double> &coords, double energy, bool accepted, MC* mc);

	virtual Array<double> get_histogram(){
		std::vector<double> vecdata =_hist->get_vecdata();
		Array<double> histogram(vecdata);
		Array<double> histogram2(histogram.copy());
		return histogram2;
	}

	virtual void print_terminal(size_t ntot){
				_hist->print_terminal(ntot);};

	virtual double get_max(){
		double max_;
		max_ = _hist->max();
		return max_;
	};

	virtual double get_min(){
			double min_;
			min_ = _hist->min();
			return min_;
		};
};

RecordEnergyHistogram::RecordEnergyHistogram(double min, double max, double bin, size_t eqsteps):
			_hist(new pele::Histogram(min, max, bin)),_bin(bin),
			_eqsteps(eqsteps),_count(0){}

void RecordEnergyHistogram::action(Array<double> &coords, double energy, bool accepted, MC* mc) {
	_count = mc->get_iterations_count();
	if (_count > _eqsteps)
		_hist->add_entry(energy);
}

}
#endif
