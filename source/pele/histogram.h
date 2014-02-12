#ifndef _HISTOGRAM_H__
#define _HISTOGRAM_H__

#include <math.h>
#include <algorithm>
#include <list>
#include "array.h"

using std::list;
using std::runtime_error;

namespace pele{

/*Dynamic histogram class that expand if energies outside of the initial bounds are found.
 * Being generous on the initial bounds saves a lot of time in reallocation of memory, at
 * the cost of memory preallocation.
 * Note:
 * floor always casts towards minus infinity
 * a list is used instead of a vector because more efficient at pushing forward
 * begin and end return list iterators to the beginning and the end of the
 * histogram respectively.
 * */

class Histogram{
protected:
	double _max, _min, _bin, _N;
	list<size_t> _hist;
public:
	Histogram(double max, double min, double bin);
	~Histogram() {}
	void add_entry(double entry);
	double max(){return _max;};
	double min(){return _min;};
	double bin(){return _bin;};
	double size(){return _N;};
	list<size_t>::iterator begin(){return _hist.begin();};
	list<size_t>::iterator end(){return _hist.end();};
};

Histogram::Histogram(double max, double min, double bin):
		_max(max),_min(min),_bin(bin),_N(floor((_max - _min) / _bin)),
		_hist(_N,0)
		{}

Histogram::add_entry(double E){
	size_t i, newlen;
	i = floor((E - _min)/_bin);
	if (i < _N && i > 0)
		_hist[i] += 1;
	else if (i > _N)
	{
		newlen = i - _N;
		for(size_t j = 0; j < newlen-1; ++j)
		{
			_hist.push_back(0);
		}
		_hist.push_back(1);
		_max = E;
		_N = floor((E - _min) / _bin);
	}
	else
	{
		newlen = i;
		for(size_t j = 0; j < newlen-1; ++j)
		{
			_hist.push_front(0);
		}
		_hist.push_front(1);
		_min = E;
		_N = floor((_max - E) / _bin);
	}
}

}

#endif
