#ifndef _PELE_LBFGS_H__
#define _PELE_LBFGS_H__

#include <vector>
#include "base_potential.h"
#include "array.h"
#include "optimizer.h"

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
	  double _finc, _fdec, _fa, _tol;
	  double _astart, _potential;
	  pele::Array<double> _coords;
	  std::vector<double> v;

    public :

	    /**
		* Constructor
		*/
	  FIRE(pele::BasePotential * potential, const pele::Array<double> & x0,
    		  double maxstep, double initdt, double dtmax,
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

  FIRE::FIRE(pele::BasePotential * potential, const pele::Array<double> & x0,
		  double maxstep, double initdt, double dtmax,
		  size_t Nmin, double finc, double fdec, double fa, double astart,
		  double iprint, double tol=1e-4):
    	  _potential(potential), _coords(x0.copy()), maxstep_(maxstep),
    	  _initdt(initdt), _dtmax(dtmax), _astart(astart), _Nmin(Nmin), _finc(finc),
    	  _fdec(fdec), _fa(fa), _tol(tol), iprint_(iprint)
  {
	_v.assign(_coords.size(),0);
	g_.assign(_coords.size(),0);
	x_.assign(_coords.size(),0);
	for(size_t i=0; i!=_coords.size(); ++i)
	{
		x_[i] = _coords[i];
	}
  }


}

#endif
