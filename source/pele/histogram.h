#ifndef _PELE_HISTOGRAM_H
#define _PELE_HISTOGRAM_H

#include <math.h>
#include <algorithm>
#include <list>
#include "array.h"

using std::vector;
using std::runtime_error;

namespace pele{

/*Dynamic histogram class that expand if energies outside of the initial bounds are found.
 * Being generous on the initial bounds saves a lot of time in reallocation of memory, at
 * the cost of memory preallocation.
 * Notes:
 * ->floor always casts towards minus infinity
 * ->a list is used instead of a vector because more efficient at pushing forward
 * ->begin and end return list iterators point to the beginning and the end of the
 *   histogram respectively.
 * */

class Histogram{
protected:
	double _max, _min, _bin;
	int _N;
	vector<size_t> _hist;
public:
	Histogram(double min, double max, double bin);
	~Histogram() {}
	void add_entry(double entry);
	double max(){return _max;};
	double min(){return _min;};
	double bin(){return _bin;};
	double size(){return _N;};
	vector<size_t>::iterator begin(){return _hist.begin();};
	vector<size_t>::iterator end(){return _hist.end();};
	void print(){
		for(int i=0; i<_hist.size();++i)
		{
			std::cout << i << "-" << (i+1) << ": ";
			std::cout << std::string(_hist[i],'*') << std::endl;
		}
	};
};

Histogram::Histogram(double min, double max, double bin):
		_max(max),_min(min),_bin(bin),_N(floor(0.5 + ((max - min) / bin))),
		_hist(_N,0)
		{
			std::cout<<"histogram is of size "<<_N<<std::endl;
		}

void Histogram::add_entry(double E){
	int i, newlen;
	i = floor(0.5+((E-_min)/_bin));
	if (i <= _N && i >= 0)
		_hist[i] += 1;
	else if (i > _N)
	{
		newlen = i - _N;
		_hist.insert(_hist.end(),(newlen-1),0);
		_hist.push_back(1);
		_max = floor((E/_bin)+0.5)*_bin; //round to nearest increment
		_N = floor(((E - _min) / _bin)+0.5);
		assert(_hist.size() == _N);
	}
	else
	{
		newlen = abs(i);
		_hist.insert(_hist.begin(),(newlen-1),0);
		_hist.insert(_hist.begin(),1);
		_min = floor((E/_bin)+0.5)*_bin; //round to nearest increment
		_N = floor(0.5+((_max-_min)/_bin));
		std::cout<<"E "<<E<<"_size "<<_hist.size()<<" i "<<i<<" bin "<<_bin<<" min "<<_min<<" max "<<_max<<"(_max - E)/_bin)= " <<floor((_max-E)/_bin)<<" newlen i "<<newlen<<std::endl;
		assert(_hist.size() == _N);
	}
}
}

#endif
