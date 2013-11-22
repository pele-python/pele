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
	  pele::Array<double> _x, _v, _f, _fold, _m, _vstart, _fstart, _xstart, _default_array;
	  pele::BasePotential * _potential
	  double _E, _dt, _Estart, _dtstart;

  	  public:

	  BaseIntegrator(pele::BasePotential * potential, pele::Array<double> x, double dt, pele::Array<double> v = _default_array,
			  pele::Array<double> m = _default_array):

				_potential(potential), _x(x), _xstart(x.copy()), _dt(dt),
	    		_dtstart(dt), _f(x.size()), _fold(x.size())

	  	  	  {
		  	  	  int i,j;
	  	  	  	  _E = _potential.energy_gradient(_x, _f);
	  	  	  	  _Estart = _E;
	  	  	  	  _fstart(_f.copy());

	  	  	  	  if (v == _default_array)
	  	  	  		  _v(x.size(),0.);
	  	  	  	  else
	  	  	  	  {
	  	  	  		  assert(v.size() == x.size());
	  	  	  		  _v(v); // NOTE: wrap velocity, it does not copy it
	  	  	  	  }

	  	  	  	  _vstart(_v.copy());

	  	  	  	  if (m == _default_array)
	  	  	  		  _m(x.size(),1.);
	  	  	  	  else
	  	  	  	  {
	  	  	  		  assert(m.size() == x.size()/3);
	  	  	  		  j = 0;
	  	  	  		  _m(m.size());

	  	  	  		  for(i=0; i != x.size(); ++i)
	  	  	  		  {
	  	  	  			  _m[i] = m[j];
	  	  	  			  if ((i+1)%3 == 0){ ++j; }
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

	  virtual void reset_dt(){ _dt = _dtstart; }

	  virtual void reset_v(){ _v(_vstart.copy()); }

	  virtual void reset_x(){ _x(_xstart.copy()); }

	  virtual void reset_f(){ _f(_fstart.copy()); }

	  virtual void reset_vfxE()
	  {
	  	_v(_vstart.copy());
	  	_x(_xstart.copy());
	  	_f(_fstart.copy());
	  	_E(_Estart.copy());
	  }

	  virtual void reset_all()
	  {
	   	_v(_vstart.copy());
	   	_x(_xstart.copy());
	   	_f(_fstart.copy());
	   	_E(_Estart.copy());
	   	_dt(_dtstart.copy());
	  }

	  virtual void wrap_v(pele::Array<double>& v){ v(_v); }

	  virtual void wrap_f(pele::Array<double>& f){ f(_f); }

	  virtual void wrap_vstart(pele::Array<double> &v){ v(_vstart); }

	  virtual void wrap_xstart(pele::Array<double> &x){ x(_xstart); }

	  virtual double get_energy(){ return _E; }
  };

}


#endif
