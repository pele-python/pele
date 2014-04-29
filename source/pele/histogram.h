#ifndef _PELE_HISTOGRAM_H
#define _PELE_HISTOGRAM_H

#include <math.h>
#include <algorithm>
#include <list>
#include "array.h"
#include <iostream>
#include <limits>

using std::vector;
using std::runtime_error;
using std::sqrt;

namespace pele{

/*Dynamic histogram class that expand if energies outside of the initial bounds are found.
 * Being generous on the initial bounds saves a lot of time in reallocation of memory, at
 * the cost of memory preallocation.
 * Notes:
 * ->floor always casts towards minus infinity
 * ->a list is used instead of a vector because more efficient at pushing forward
 * ->begin and end return list iterators point to the beginning and the end of the
 *   histogram respectively.
 * -> the most basic test that histogram must satisfy is that there must be as many
 * 	 beads as the number of iterations (commented out at the end of the script)
 * */

class Histogram{
protected:
	double _max, _min, _bin, _eps;
	int _N;
	vector<double> _hist;
public:
	int _niter;
	Histogram(double min, double max, double bin);
	~Histogram() {}
	void add_entry(double entry);
	double max(){return _max;};
	double min(){return _min;};
	double bin(){return _bin;};
	size_t size(){return _N;};
	vector<double>::iterator begin();
	vector<double>::iterator end();
	vector<double> get_vecdata(){return _hist;};
	void print_terminal(size_t ntot){
		for(size_t i=0; i<_hist.size();++i)
		{
			std::cout << i << "-" << (i+1) << ": ";
			std::cout << std::string(_hist[i]*10000/ntot,'*') << std::endl;
		}
	};
	void resize(double E, int i);
};

Histogram::Histogram(double min, double max, double bin):
		_max(floor((max/bin)+1)*bin),_min(floor((min/bin))*bin),_bin(bin),
		_eps(std::numeric_limits<double>::epsilon()),_N((_max - _min) / bin),
		_hist(_N,0),_niter(0)
		{
        #ifdef DEBUG
			std::cout<<"histogram is of size "<<_N<<std::endl;
        #endif
		}

inline void Histogram::add_entry(double E){
	int i;
	E = E + _eps; //this is a dirty hack, not entirely sure of its generality and possible consequences, tests seem to be fine
	i = floor((E-_min)/_bin);
	if (i < _N && i >= 0)
	{
		_hist[i] += 1;
		++_niter;
	}
	else
		this->resize(E,i);

	/*THIS IS A TEST*/
	/*int renorm = 0;
	 * for(vector<size_t>::iterator it = _hist.begin();it != _hist.end();++it)
	  {
		  renorm += *it;
	  }

	if (renorm != _niter)
	{
		std::cout<<" E "<<E<<"\n niter "<<_niter<<"\n renorm "<<renorm<<"\n min "<<_min<<"\n max "<<_max<<"\n i "<<i<<"\n N "<<_N<<std::endl;
		assert(renorm == _niter);
	}*/
}

inline void Histogram::resize(double E, int i){
	int newlen;
	if (i >= _N)
		{
			newlen = (i + 1) - _N;
			_hist.insert(_hist.end(),(newlen-1),0);
			_hist.push_back(1);
			++_niter;
			_max = floor((E/_bin)+1)*_bin; //round to nearest increment
			_N = round((_max - _min) / _bin); //was round
			if ( (int) _hist.size() != _N)
			{
				std::cout<<" E "<<E<<"\n niter "<<_niter<<"\n size "<<_hist.size()<<"\n min "<<_min<<"\n max "<<_max<<"\n i "<<i<<"\n N "<<_N<<std::endl;
				assert( (int) _hist.size() == _N);
				exit (EXIT_FAILURE);
			}
			std::cout<<"resized above at niter "<<_niter<<std::endl;
		}
	else if (i < 0)
		{
			newlen = -1*i;
			_hist.insert(_hist.begin(),(newlen-1),0);
			_hist.insert(_hist.begin(),1);
			++_niter;
			_min = floor((E/_bin))*_bin; //round to nearest increment
			_N = round((_max-_min)/_bin); //was round
			if ( (int) _hist.size() != _N)
			{
				std::cout<<" E "<<E<<"\n niter "<<_niter<<"\n size "<<_hist.size()<<"\n min "<<_min<<"\n max "<<_max<<"\n i "<<i<<"\n N "<<_N<<std::endl;
				assert( (int) _hist.size() == _N);
				exit (EXIT_FAILURE);
			}
			std::cout<<"resized below at niter "<<_niter<<std::endl;
		}
	else
		{
			std::cerr<<"histogram encountered unexpected condition"<<std::endl;
			std::cout<<" E "<<E<<"\n niter "<<_niter<<"\n min "<<_min<<"\n max "<<_max<<"\n i "<<i<<"\n N "<<_N<<std::endl;
		}
}

vector<double>::iterator Histogram::begin(){return _hist.begin();}
vector<double>::iterator Histogram::end(){return _hist.end();}

}

#endif
