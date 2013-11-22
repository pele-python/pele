#ifndef _PELE_INTEGRATOR_H__
#define _PELE_INTEGRATOR_H__

#include <assert.h>
#include <vector>
#include <stdexcept>
#include <iostream>
#include "array.h"

namespace pele{

class BaseIntegrator
  {
  	  protected:
	  std::vectror<double>& _x;
	  std::vector<double> _v, _f, _fold, _m, _vstart, _fstart, _xstart, _default_vector;
	  pele::BasePotential * _potential
	  double _E, _dt, _Estart, _dtstart;

  	  public:

	  BaseIntegrator(pele::BasePotential * potential, std::vector<double>& x, double dt, std::vector<double> m = _default_vector):
	    		_potential(potential), _x(x), _xstart(x), _dt(dt),
	    		_dtstart(dt), _f(x.size(),0), _fold(x.size(),0),
	    		_v(x.size(),0), _vstart(x.size(),0)
	    		{
	  	  	  	  _E = _potential.energy_gradient(_x, _f);
	  	  	  	  _fstart(_f);

	  	  	  	  if (m == _default_vector)
	  	  	  		  _m(x.size(), 1.);
	  	  	  	  else
	  	  	  	  {
	  	  	  		  assert(m.size() == x.size()/3);
	  	  	  		  int j = 0;
	  	  	  		  _m(x);

	  	  	  		  for(int i=0; i != x.size(); ++i)
	  	  	  		  {
	  	  	  			  if (i%3 == 0){ ++j; }
	  	  	  			  _m[i] = m[j];
	  	  	  		  }
	  	  	  	  }
	    		}

	  virtual ~BaseIntegrator() {}

	  virtual void oneiteration()
	  {
		  throw std::runtime_error("BaseIterator::oneiteration must be overloaded");
	  }

	  virtual void run(int N)
	  {
		  throw std::runtime_error("BaseIterator::run must be overloaded");
	  }

	  /*reset initial velocities*/

	  virtual void set_dt(double newdt){ _dt = newdt; }

	  virtual void reset_dt(){ _dt(_dtstart); }

	  virtual void reset_v(){ _v(_vstart); }

	  virtual void reset_x(){ _x(_xstart); }

	  virtual void reset_f(){ _f(_fstart); }

	  virtual void reset_vfxE()
	  {
	  	_v(_vstart);
	  	_x(_xstart);
	  	_f(_fstart);
	  	_E(_Estart);
	  }

	  virtual void get_v(std::vector<double> &v){ v = _v; }

	  virtual void get_f(std::vector<double> &f){ f = _f; }

	  virtual double get_energy(){ return _E; }
  };

}


#endif
