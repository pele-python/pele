#ifndef _PELE_VELOCITY_VERLET_H__
#define _PELE_VELOCITY_VERLET_H__

#include "base_integrator.h"

namespace pele {

	class VelocityVerlet: public BaseIntegrator
	  {
	  private:
		  std::vector<double> _default_vector;
	  public:
		  /*Constructor*/
		  VelocityVerlet(pele::BasePotential * potential, std::vector<double>& x, double dt, std::vector<double> m = _default_vector){}
	  };

	  VelocityVerlet::VelocityVerlet(pele::BasePotential * potential, std::vector<double> x, double dt, std::vector<double> m = _default_vector):
			  BaseIntegrator(potential, x, dt, m)
		{}

	  void VelocityVerlet::oneiteration()
	  {
		  int j = 0;

		  for(int i =0; i != _x.size(); ++i)
		  {
			_x[i] += _dt * (_v[i] + 0.5 * _dt * _f[i] / _m[i]);

			_fold(_f);

			_E = _potential.energy_gradient(_x, _f);

			_v[i] += 0.5 * _dt * (_fold[i] + _f[i]) / _m[i];

		  }
	  }

	  void VelocityVerlet::run(int N)
	  {
		 for(int i =0; i != N; ++i)
		 {
			 this.oneiteration();
		 }
	  }

}
