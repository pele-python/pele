#ifndef _PELE_ACTIONS_H
#define _PELE_ACTIONS_H

#include <math.h>
#include <algorithm>
#include <list>
#include "array.h"
#include "mc.h"
#include "histogram.h"

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
public:
	RecordEnergyHistogram(double min, double max, double bin);
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

RecordEnergyHistogram::RecordEnergyHistogram(double min, double max, double bin):
			_hist(new pele::Histogram(min, max, bin)),_bin(bin){}

void RecordEnergyHistogram::action(Array<double> &coords, double energy, bool accepted, MC* mc) {
		_hist->add_entry(energy);
}

/*Record displacement square histogram
*/

class RecordDisp2Histogram : public RecordEnergyHistogram {
protected:
	pele::Array<double> _origin, _rattlers, _distance;
	size_t _N;
public:
	RecordDisp2Histogram(pele::Array<double> origin, pele::Array<double> rattlers, double min, double max, double bin);
	virtual ~RecordDisp2Histogram() {delete _hist;}
	virtual void action(Array<double> &coords, double energy, bool accepted, MC* mc);
};

RecordDisp2Histogram::RecordDisp2Histogram(pele::Array<double> origin, pele::Array<double> rattlers, double min, double max, double bin):
		RecordEnergyHistogram(min, max, bin),
		_origin(origin.copy()),_rattlers(rattlers.copy()),
		_distance(origin.size()),_N(origin.size()){}

void RecordDisp2Histogram::action(Array<double> &coords, double energy, bool accepted, MC* mc) {
		double _d;

		_distance.assign(coords);
		//compute distances subtracting the origin's coordinates
		_distance -= _origin;
		//set to 0 distances of rattlers
		for (size_t j = 0; j < _N; ++j){
			_distance[j] *= _rattlers[j];
		}
		//compute square displacement from origin
		_d = norm(_distance);
		_hist->add_entry(_d*_d);
}









}
#endif
