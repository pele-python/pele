#ifndef _PELE_MODIFIED_FIRE_H__
#define _PELE_MODIFIED_FIRE_H__

#include <vector>
#include "base_potential.h"
#include "array.h"
#include "optimizer.h"
#include "base_integrator.h"

using std::vector;

namespace pele{
  /**
   * An implementation of the *modified* FIRE optimization algorithm in c++.
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

	class MODIFIED_FIRE : public GradientOptimizer{
    private :
	  double _Nmin, _dt, _dtmax;
	  double _finc, _fdec, _fa;
	  double _astart, _a;
	  pele::BasePotential * _potential;
	  pele::BaseIntegrator * _integrator;
	  pele::Array<double> _xstart, _f, _v, _x, _xold;

    public :

	    /**
		* Constructor
		*/
	  MODIFIED_FIRE(pele::BasePotential * potential, pele::BaseIntegrator * integrator,
			  const pele::Array<double> & x0, double maxstep, double dtstart, double dtmax,
    		  size_t Nmin, double finc, double fdec, double fa, double astart,
    		  double iprint, double tol=1e-4) {}
      /**
       * Destructorgit undo rebase
       */
      ~MODIFIED_FIRE() {}

      /**
       * Do one iteration iteration of the optimization algorithm
       */

      void one_iteration();

  };

  MODIFIED_FIRE::MODIFIEDFIRE(pele::BasePotential * potential, pele::BaseIntegrator * integrator,
		  const pele::Array<double> & x0, double maxstep, double dtstart, double dtmax,
		  size_t Nmin, double finc, double fdec, double fa, double astart,
		  double iprint, double tol=1e-4):
		  GradientOptimizer(potential, x0, tol),							//call to GradienOptimizer constructor to initialise inherited variables
		  _integrator(integrator(potential, x0, dtstart)),						//call to BaseIntegrator constructor to initialise integrator variables
		  _xstart(x0.copy()), _dtstart(dtstart),
		  _dtmax(dtmax), _astart(astart), _a(astart),
		  _Nmin(Nmin), _finc(finc), _fdec(fdec),
		  _fa(fa), _f(x0.size()), _v(x0.size()),
    	  _potential(potential), _dt(dtstart),
    	  _xold(x0.copy())

  {
	  set_maxstep(maxstep);
	  set_iprint(iprint);
	  _integrator.wrapv(_v); 	//the velocity array wraps the integrator velocity array so that it updates concurrently
	  _integrator.wrapf(_f); 	//the force array wraps the integrator force array so that it updates concurrently
	  _integrator.wrapx(_x);	//the coordinates array wraps the integrator coordinates array so that it updates concurrently
  }

  MODIFIED_FIRE::one_iteration()
  {
	  _integrator.oneiteration();
	  double P = arraydot(_v,_f);

	  if (P > 0)
	  {
		  /*equation written in this conditional statement _v = (1- _a)*_v + _a * funit * vnorm*/

		  /*double vnorm = arraynorm(_v);
		  pele::array<double> funit(_f.copy());
		  arrayunit(funit); 					//transforms _f in a unit vector
		  arrayscalprod(funit, _a*vnorm); 		//_funit -> _a * funit * vnorm
		  arrayscalprod(_v, 1.-_a);				// _v -> (1-_a) * _v
		  arraysum(_v,funit,_v);				// _v = (1- _a)*_v + _a * funit * vnorm*/

		  double ifnorm = 1. / arraynorm(_f);
		  double ivnorm = 1. / arraynorm(_v);

		  for (size_t i =0; i < v.size(); ++v)
		  {
			  v[i] += (1. - _a) * _v[i] + _a * _f[i] * _v[i] * ifnorm * ivnorm;
		  }

		  if (iter_number_ > _Nmin)
		  {
			  _dt = std::min(_dt* _finc, _dtmax);
			  integrator.set_dt(_dt);
			  _a *= _fa;
		  }
		  iter_number_ += 1;
		  _xold(_x.copy());
	  }
	  else
	  {
		  _dt *= _fdec;
		  integrator.set_dt(_dt);
		  integrator.set_v_to_vstart(); 		//reset velocity to initial (0)
		  _x(_xold.copy());						//reset position to the one before the step (core of modified fire)
		  _a = _astart;
		  iter_number_ = 0;
	  }
  }
}

#endif
