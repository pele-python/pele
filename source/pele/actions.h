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

/*Adjust Step*/

class AdjustStep : public Action {
protected:
	double _factor;
public:
	AdjustStep(double factor);
	virtual ~AdjustStep() {}
	virtual void action(Array<double> &coords, double energy, bool accepted, MC* mc);
};

AdjustStep::AdjustStep(double factor):
			_factor(factor){}


void AdjustStep::action(Array<double> &coords, double energy, bool accepted, MC* mc) {
		if (accepted == false)
			mc->_stepsize *= _factor;
		else
			mc->_stepsize /= _factor;
	}

/*Record energy histogram
 * Note:
 * the histogram has to be initialised elsewhere and it is passed by reference
 * */

class RecordEnergyHistogram : public Action {
protected:
	pele::Histogram * _hist;
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
};

RecordEnergyHistogram::RecordEnergyHistogram(double min, double max, double bin):
			_hist(new pele::Histogram(min, max, bin)){}

void RecordEnergyHistogram::action(Array<double> &coords, double energy, bool accepted, MC* mc) {
		_hist->add_entry(energy);
}



}
#endif
