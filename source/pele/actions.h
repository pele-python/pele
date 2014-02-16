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
	pele::Histogram* _hist;
public:
	RecordEnergyHistogram(pele::Histogram * hist);
	virtual ~RecordEnergyHistogram() {delete _hist;}
	virtual void action(Array<double> &coords, double energy, bool accepted, MC* mc);
};

RecordEnergyHistogram::RecordEnergyHistogram(pele::Histogram * hist):
			_hist(hist){}

void RecordEnergyHistogram::action(Array<double> &coords, double energy, bool accepted, MC* mc) {
		_hist->add_entry(energy);
}



}
#endif
