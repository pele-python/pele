#ifndef _PELE_VELOCITY_VERLET_H__
#define _PELE_VELOCITY_VERLET_H__

#include "base_integrator.h"

namespace pele {

	class VelocityVerlet: public BaseIntegrator
	  {
	  public:
		  /*Constructor*/
		VelocityVerlet(pele::BasePotential * potential, pele::Array<double> x, double dt,
					pele::Array<double> v = pele::Array<double>(),pele::Array<double> g = pele::Array<double>(),
					pele::Array<double> m = pele::Array<double>());

		void oneiteration()
		  	  {
		  		   /* the minuses in the following expressions are due to the fact that
		  		   * the gradients rather than the forces appear in the expression
		  		   */
			  	  size_t i;

			  	  for(i =0; i < _x.size(); ++i)
		  		  {

					_gold[i] = _g[i]; //save gradient as old g

					_x[i] += _dt * (_v[i] - 0.5 * _dt * _g[i] / _m[i]);	//update position

					*_E = (*_potential).get_energy_gradient(_x, _g);	//update gradient

					_v[i] -= 0.5 * _dt * (_gold[i] + _g[i]) / _m[i]; 	//update velocity

		  		  }
		  	  }

		void run(size_t const N)
		  	  {
		  		 for(size_t i =0; i < N; ++i)
		  		 {
		  			 oneiteration();
		  		 }
		  	  }

	  };

	VelocityVerlet::VelocityVerlet(pele::BasePotential * potential, pele::Array<double> x, double dt,
			pele::Array<double> v, pele::Array<double> g, pele::Array<double> m):
			BaseIntegrator(potential, x, dt, v, g, m) //initialise base integrator from which this class is inherited
			{}

}
#endif
