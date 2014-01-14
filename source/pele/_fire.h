#ifndef _PELE_LBFGS_H__
#define _PELE_LBFGS_H__

#include <vector>
#include "base_potential.h"
#include "array.h"
#include "optimizer.h"
#include "base_integrator.h"

using std::vector;

namespace pele{
  /**
   * An implementation of the FIRE optimization algorithm in c++.
   *
   * The Fast Inertial Relaxation Engine is an optimization algorithm based
   * on molecular dynamics with modifications to the velocity and adaptive
   * time steps. The method is based on a blind skier searching for the
   * bottom of a valley and is described and tested here:
   *
   * Erik Bitzek, Pekka Koskinen, Franz Gaehler, Michael Moseler, and Peter Gumbsch.
   * Phys. Rev. Lett. 97, 170201 (2006)
   * http://link.aps.org/doi/10.1103/PhysRevLett.97.170201
   */

	class FIRE : public GradientOptimizer{
    private :
	  double _Nmin, _dt, _dtmax;
	  double _finc, _fdec, _fa;
	  double _astart;
	  pele::BasePotential * _potential;
	  pele::BaseIntegrator * _integrator;
	  pele::Array<double> _xstart, _f, _v;

    public :

	    /**
		* Constructor
		*/
	  FIRE(pele::BasePotential * potential, pele::BaseIntegrator * integrator,
			  const pele::Array<double> & x0, double maxstep, double dtstart, double dtmax,
    		  size_t Nmin, double finc, double fdec, double fa, double astart,
    		  double iprint, double tol=1e-4) {}
      /**
       * Destructor
       */
      ~FIRE() {}

      /**
       * Do one iteration iteration of the optimization algorithm
       */

      void one_iteration();

  };

  FIRE::FIRE(pele::BasePotential * potential, pele::BaseIntegrator * integrator,
		  const pele::Array<double> & x0, double maxstep, double dtstart, double dtmax,
		  size_t Nmin, double finc, double fdec, double fa, double astart,
		  double iprint, double tol=1e-4):
		  GradientOptimizer(potential, x0, tol=1e-4),
		  _xstart(x0.copy()), _dtstart(dtstart), _dtmax(dtmax),
		  _astart(astart), _Nmin(Nmin), _finc(finc),
    	  _fdec(fdec), _fa(fa), _f(x0.size()), _v(x0.size())

  {
	  set_maxstep(maxstep);
	  set_iprint(iprint);
  }

}

#endif
