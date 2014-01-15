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
		/* the following are the arrays:
		 *
		 * _x: position
		 * _v: velocity
		 * _g: gradient
		 * _gold: gradient at previous time step
		 * _m: masses
		 * _v: initial velocities
		 * _gstart: initial gradient
		 * _xstart: initial positions
		 * _default_array: an empty array
		 *
		 *these are the constants:
		 *
		 *_E: energy
		 *_dt: time step
		 *_Estart: initial energy
		 *_dtstart: initial time step
		 */
	  pele::Array<double> _x, _v, _g, _gold, _m, _vstart, _gstart, _xstart, _default_array;
	  pele::BasePotential * _potential
	  double _dt, _Estart, _dtstart;
	  double* _E;

  	  public:

	  	 /*BaseIntegrator constructor, assigns value to the protected arrays and constants*/
	  BaseIntegrator(pele::BasePotential * potential, pele::Array<double>& x, double dt, pele::Array<double>& v = _default_array,
			  pele::Array<double>& g = _default_array, pele::Array<double>& m = _default_array):

				_potential(potential), _x(x), _xstart(x.copy()), _dt(dt),
	    		_dtstart(dt), _gold(x.size())

	  	  	  {
		  	  	  int i,j;

		  	  	  if (g.empty())
		  	  		  _g(x.size());
		  	  	  else
		  	  	  {
		  	  		assert(g.size() == x.size());
		  	  		_g(g); // NOTE: wrap gradient, it does not copy it
		  	  	  }

		  	  	  *_E = _potential.energy_gradient(_x, _g); //potential.energy_gradient returns the energy and modifies the gradient vector by reference
	  	  	  	  _Estart = *_E;
	  	  	  	  _gstart(_g.copy());

	  	  	  	  if (v.empty())
	  	  	  		  _v(x.size(),0.);
	  	  	  	  else
	  	  	  	  {
	  	  	  		  assert(v.size() == x.size());
	  	  	  		  _v(v); // NOTE: wrap velocity, it does not copy it
	  	  	  	  }

	  	  	  	  _vstart(_v.copy());

	  	  	  	  if (m.empty())
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

	  virtual void reset_dt(){ _dt = _dtstart; }

	  virtual void reset_v(){ _v(_vstart.copy()); }

	  virtual void reset_x(){ _x(_xstart.copy()); }

	  virtual void reset_g(){ _g(_gstart.copy()); }

	  virtual void wrapx(pele::Array<double>& x){ x(_x); }

	  virtual void wrapv(pele::Array<double>& v){ v(_v); }

	  virtual void wrapg(pele::Array<double>& g){ g(_g); }

	  virtual void wrapE(double E){ _E = &E;}

	  virtual void wrap_vstart(pele::Array<double> &v){ v(_vstart); }

	  virtual void wrap_xstart(pele::Array<double> &x){ x(_xstart); }

	  virtual double get_energy(){ return *_E; }
  };

}


#endif
