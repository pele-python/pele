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
				  pele::Array<double>& f = _default_array, pele::Array<double>& m = _default_array){}
	  };

	  VelocityVerlet::VelocityVerlet(pele::BasePotential * potential, pele::Array<double> x, double dt, pele::Array<double>& v = _default_array,
			  pele::Array<double>& f = _default_array, pele::Array<double>& m = _default_array):
			  BaseIntegrator(potential, x, dt, v, f, m)
		{}

	  void VelocityVerlet::oneiteration()
	  {
		  int j = 0;

		  for(int i =0; i < _x.size(); ++i)
		  {
			_x[i] += _dt * (_v[i] + 0.5 * _dt * _f[i] / _m[i]);

			_fold(_f);

			_E = _potential.energy_gradient(_x, _f);

			_v[i] += 0.5 * _dt * (_fold[i] + _f[i]) / _m[i];

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
