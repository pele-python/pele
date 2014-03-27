#ifndef _PELE_INTEGRATOR_H__
#define _PELE_INTEGRATOR_H__

#include <assert.h>
#include <vector>
#include <stdexcept>
#include <iostream>
#include "array.h"
#include "base_potential.h"

namespace pele{

/*base class for numerical integration methods for the solution of ordinary differential equations,
 * these could be Euler method, Verlet, Velocity verlet, etc. */

class BaseIntegrator
  {
  	  protected:
		/* the following are the arrays:
		 *
		 * _x: position
		 * _v: velocity
		 * _g: gradient
		 * _gold: gradient at previous time step
		 * _dx: displacement vector (needed to control step size)
		 * _m: masses
		 * _v: initial velocities
		 *
		 *these are the variables:
		 *
		 *_E: energy
		 *_dt: time step
		 *_maxstep: maximum step to take
		 *_Estart: initial energy
		 *_dtstart: initial time step
		 */

	  pele::BasePotential * _potential;
	  pele::Array<double> _x;
	  double _dt, _dtstart, _Estart, _maxstep;
	  pele::Array<double> _gold, _v, _g, _m, _dx;
	  double * _E;

  	  public:

	  /*BaseIntegrator constructor, assigns value to the protected arrays and constants*/
	  BaseIntegrator(pele::BasePotential * potential, pele::Array<double> x, double dt, double maxstep, pele::Array<double> &v,
			  pele::Array<double> &g, pele::Array<double> &m);

	  virtual ~BaseIntegrator() {}

	  virtual inline void oneiteration()
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
		  if (_x.size() != x.size())
			  std::cout<<"warning: wrapping positions arrays of different sizes, this is dangerous.\n";
		  if (x.reference_count() > 1)
			  std::cout<<"warning: x reference count > 1, this is dangerous.\n";
		  x.wrap(_x);
	  }

	  void wrap_v(pele::Array<double>& v){
		  if(_v.size() != v.size())
			  std::cout<<"warning: wrapping velocity arrays of different sizes, this is dangerous.\n";
		  if (v.reference_count() > 1)
			  std::cout<<"warning: v reference count > 1, this is dangerous.\n";
		  v.wrap(_v);
	  }

	  void wrap_g(pele::Array<double>& g){
		  if(_g.size() != g.size())
			  std::cout<<"warning: wrapping gradient arrays of different sizes, this is dangerous.\n";
		  if (g.reference_count() > 1)
			  std::cout<<"warning: g reference count > 1, this is dangerous.\n";
		  g.wrap(_g);
	  }

	  void wrap_gold(pele::Array<double>& gold){
		  if(_gold.size() != gold.size())
			  std::cout<<"warning: wrapping gold gradient arrays of different sizes, this is dangerous.\n";
		  if (gold.reference_count() > 1)
			  std::cout<<"warning: gold reference count > 1, this is dangerous.\n";
		  gold.wrap(_gold);
	  }

	  void wrap_E(double &E){
		  double Eval = *_E;
		  _E = &E; //relocate address
		  *_E = Eval;
	  }

	  double get_energy(){ return *_E; }
  };

	BaseIntegrator::BaseIntegrator(pele::BasePotential * potential, pele::Array<double> x, double dt, double maxstep,
			pele::Array<double> &v, pele::Array<double> &g, pele::Array<double> &m):
				_potential(potential), _x(x), _dt(dt), _dtstart(dt),
				_maxstep(maxstep), _gold(x.size()), _dx(x.size()),
				_E(new double)

	  	  	  {
		  	  	  size_t i,j;
		  	  	  size_t k;

		  	  	  if (g.empty())
		  	  	  {
		  	  		  _g.resize(x.size());
		  	  	  }
		  	  	  else
		  	  	  {
		  	  		assert(g.size() == x.size());
		  	  		_g.wrap(g); // NOTE: wrap gradient, it does not copy it
		  	  	  }

		  	  	  *_E = _potential->get_energy_gradient(_x, _g); //potential.energy_gradient returns the energy and modifies the gradient vector by reference
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
