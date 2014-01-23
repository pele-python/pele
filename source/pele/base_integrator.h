#ifndef _PELE_INTEGRATOR_H__
#define _PELE_INTEGRATOR_H__

#include <assert.h>
#include <vector>
#include <stdexcept>
#include <iostream>
#include "array.h"
#include "base_potential.h"

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
		 * _default_array: an empty array
		 *
		 *these are the constants:
		 *
		 *_E: energy
		 *_dt: time step
		 *_Estart: initial energy
		 *_dtstart: initial time step
		 */

	  pele::BasePotential * _potential;
	  pele::Array<double> _x;
	  double _dt, _dtstart, _Estart;
	  pele::Array<double> _gold, _v, _g, _m;
	  double* _E;

  	  public:

	  /*BaseIntegrator constructor, assigns value to the protected arrays and constants*/
	  BaseIntegrator(pele::BasePotential * potential, pele::Array<double>& x, double dt, pele::Array<double> v = pele::Array<double>(),
			  pele::Array<double> g = pele::Array<double>(), pele::Array<double> m = pele::Array<double>());

	  virtual ~BaseIntegrator() {}

	  virtual void oneiteration()
	  {
		  throw std::runtime_error("BaseIterator::oneiteration must be overloaded");
	  }

	  virtual void run(size_t const N)
	  {
		  throw std::runtime_error("BaseIterator::run must be overloaded");
	  }

	  /*reset initial velocities*/

	  void set_dt(double newdt){ _dt = newdt; }

	  void reset_dt(){ _dt = _dtstart; }

	  void wrap_x(pele::Array<double>& x){
		  assert(_x.size() == x.size());
		  _x.wrap(x);
	  }

	  void wrap_v(pele::Array<double>& v){
		  assert(_v.size() == v.size());
		  _v.wrap(v);
	  }

	  void wrap_g(pele::Array<double>& g){
		  assert(_g.size() == g.size());
		  _g.wrap(g);
	  }

	  void wrap_gold(pele::Array<double>& gold){
		  assert(_gold.size() == gold.size());
		  _gold.wrap(gold);
	  }

	  void wrap_E(double E){ _E = &E;}

	  double get_energy(){ return *_E; }
  };

	BaseIntegrator::BaseIntegrator(pele::BasePotential * potential, pele::Array<double>& x, double dt,
			pele::Array<double> v, pele::Array<double> g, pele::Array<double> m):
				_potential(potential), _x(x), _dt(dt),
	    		_dtstart(dt), _gold(x.size())

	  	  	  {
		  	  	  size_t i,j;
		  	  	  size_t k;

		  	  	  if (g.empty())
		  	  		  _g.resize(x.size());
		  	  	  else
		  	  	  {
		  	  		assert(g.size() == x.size());
		  	  		_g.wrap(g); // NOTE: wrap gradient, it does not copy it
		  	  	  }

		  	  	  *_E = (*_potential).get_energy_gradient(_x, _g); //potential.energy_gradient returns the energy and modifies the gradient vector by reference
	  	  	  	  _Estart = *_E;

	  	  	  	  for(k=0; k<_g.size();++k)
	  	  	  	  	  {
	  	  	  	  	  	  _gold[k] = _g[k];
	  	  	  	  	  }

	  	  	  	  if (v.empty())
	  	  	  	  {
	  	  	  		  _v.resize(x.size());
	  	  	  		  for(k=0; k<x.size();++k)
	  	  	  		  	  {
	  	  	  			  	  _v[k] = 0;
	  	  	  		  	  }
	  	  	  	  }
	  	  	  	  else
	  	  	  	  {
	  	  	  		  assert(v.size() == x.size());
	  	  	  		  _v.wrap(v); // NOTE: wrap velocity, it does not copy it
	  	  	  	  }

	  	  	  	  if (m.empty())
	  	  	  	  {
	  	  	  		  _m.resize(x.size());
	  	  	  	  	  for(k=0; k<x.size();++k)
	  	  	  	  	  	  {
	  	  	  	  	  	  	  _m[k] = 1.;
	  	  	  	  	  	  }
	  	  	  	  }
	  	  	  	  else
	  	  	  	  {
	  	  	  		  assert(m.size() == x.size()/3);
	  	  	  		  j = 0;
	  	  	  		  _m.resize(m.size());

	  	  	  		  for(i=0; i != x.size(); ++i)
	  	  	  		  {
	  	  	  			  _m[i] = m[j];
	  	  	  			  if ((i+1)%3 == 0){ ++j; }
	  	  	  		  }
	  	  	  	  }
	    		}

}


#endif
