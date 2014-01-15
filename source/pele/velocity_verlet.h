#ifndef _PELE_VELOCITY_VERLET_H__
#define _PELE_VELOCITY_VERLET_H__

#include "base_integrator.h"

namespace pele {

	class VelocityVerlet: public BaseIntegrator
	  {
	  private:
		  pele::Array<double> _default_array;
	  public:
		  /*Constructor*/
		  VelocityVerlet(pele::BasePotential * potential, pele::Array<double> x, double dt, pele::Array<double>& v = _default_array,
				  pele::Array<double>& g = _default_array, pele::Array<double>& m = _default_array){}
	  };

	  VelocityVerlet::VelocityVerlet(pele::BasePotential * potential, pele::Array<double> x, double dt, pele::Array<double>& v = _default_array,
			  pele::Array<double>& g = _default_array, pele::Array<double>& m = _default_array):
			  BaseIntegrator(potential, x, dt, v, g, m) //initialise base integrator from which this class is inherited
		{}

	  void VelocityVerlet::oneiteration()
	  {
		  int j = 0;

		  /* the minuses in the following expressions are due to the fact that
		   * the gradients rather than the forces appear in the expression
		   */

		  for(int i =0; i < _x.size(); ++i)
		  {
			_x[i] -= _dt * (_v[i] + 0.5 * _dt * _g[i] / _m[i]);	//update position

			_gold(_g.copy());

			_E = _potential.energy_gradient(_x, _g); 			//update gradient

			_v[i] -= 0.5 * _dt * (_gold[i] + _g[i]) / _m[i]; 	//update velocity

		  }
	  }

	  void VelocityVerlet::run(int const N)
	  {
		 for(int i =0; i < N; ++i)
		 {
			 this.oneiteration();
		 }
	  }

}
