#ifndef _PELE_VELOCITY_VERLET_H__
#define _PELE_VELOCITY_VERLET_H__

#include "base_integrator.h"

namespace pele {

	class VelocityVerlet: public BaseIntegrator
	  {
	  public:
		  /*Constructor*/
		VelocityVerlet(pele::BasePotential * potential, pele::Array<double> x, double dt, double maxstep,
					pele::Array<double> v = pele::Array<double>(), pele::Array<double> g = pele::Array<double>(),
					pele::Array<double> m = pele::Array<double>());

		void oneiteration()
		  	  {
				   /* the minuses in the following expressions are due to the fact that
				   * the gradients rather than the forces appear in the expression
				   */
				  size_t i;
				  double normdx;

				for(i =0; i < _x.size(); ++i)
				  {
					_gold[i] = _g[i]; //save gradient as old g
					_v[i] -= 0.5 * _dt * (_gold[i] + _g[i]) / _m[i]; 	//update velocity
					_dx[i] = _dt * (_v[i] - 0.5 * _dt * _g[i] / _m[i]);	//build displacement vector
				  }

				  normdx = norm(_dx);

				  if(normdx > _maxstep){
					_dx *= (_maxstep / normdx); //resize displacement vector is greater than _maxstep
					}

				  _x += _dx;

				  *_E = _potential->get_energy_gradient(_x, _g);	//update gradient

		  	  }

		void run(size_t const N)
		  	  {
		  		 for(size_t i =0; i < N; ++i)
		  		 {
		  			 oneiteration();
		  		 }
		  	  }

	  };

	VelocityVerlet::VelocityVerlet(pele::BasePotential * potential, pele::Array<double> x, double dt, double maxstep,
			pele::Array<double> v, pele::Array<double> g, pele::Array<double> m):
			BaseIntegrator(potential, x, dt, maxstep, v, g, m) //initialise base integrator from which this class is inherited
			{}

}
#endif
