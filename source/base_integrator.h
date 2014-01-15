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

	  BaseIntegrator(pele::BasePotential * potential, pele::Array<double>& x, double dt, pele::Array<double>& v = _default_array,
			  pele::Array<double>& f = _default_array, pele::Array<double>& m = _default_array):

				_potential(potential), _x(x), _xstart(x.copy()), _dt(dt),
	    		_dtstart(dt), _fold(x.size())

	  	  	  {
		  	  	  int i,j;

		  	  	  if (f == _default_array)
		  	  		  _f(x.size());
		  	  	  else
		  	  	  {
		  	  		assert(f.size() == x.size());
		  	  		_f(f); // NOTE: wrap force, it does not copy it
		  	  	  }

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

	  virtual void run(int const N)
	  {
		  throw std::runtime_error("BaseIterator::run must be overloaded");
	  }

	  /*reset initial velocities*/

	  virtual void set_dt(double newdt){ _dt = newdt; }

	  virtual void set_x(pele::Array<double> &x){ _x(x); }

	  virtual void set_v(pele::Array<double> &v){ _v(v); }

	  virtual void set_f(pele::Array<double> &f){ _f(f); }

	  virtual void set_dtstart_to_t(){ _dtstart(_dt.copy()); }

	  virtual void set_xstart_to_x(){ _xstart(_x.copy()); }

	  virtual void set_vstart_to_v(){ _vstart(_v.copy()); }

	  virtual void set_fstart_to_f(){ _fstart(_f.copy()); }

	  virtual void set_dt_to_dtstart(){ _dt = _dtstart; }

	  virtual void set_v_to_vstart(){ _v(_vstart.copy()); }

	  virtual void set_x_to_xstart(){ _x(_xstart.copy()); }

	  virtual void set_f_to_fstart(){ _f(_fstart.copy()); }

	  virtual void set_vfxE_to_start()
	  {
	  	_v(_vstart.copy());
	  	_x(_xstart.copy());
	  	_f(_fstart.copy());
	  	_E(_Estart.copy());
	  }

	  virtual void reset_all_to_start()
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

	  virtual double get_dt(){ return _dt; }
  };

}


#endif
