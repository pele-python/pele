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

	class MODIFIED_FIRE : public Optimizer{
    private :
	  double _Nmin, _dt, _dtmax;
	  double _finc, _fdec, _fa;
	  double _astart, _a;
	  pele::BaseIntegrator * _integrator;
      pele::Array<double> _v, _xold;
      int _fire_iter_number = 0;
      int _N = 0;

    public :

	    /**
		* Constructor
		*/
	  MODIFIED_FIRE(pele::BasePotential * potential, pele::BaseIntegrator * integrator,
			  const pele::Array<double> & x0, double dtstart, double dtmax,
    		  size_t Nmin, double finc, double fdec, double fa, double astart, double tol=1e-2) {}
      /**
       * Destructorgit undo rebase
       */
      ~MODIFIED_FIRE() {}

      /**
       * Do one iteration iteration of the optimization algorithm
       */

      void one_iteration();

      /**
       * Overload initialize_func_gradient from parent class (here we need to wrap the integrator)
       */

      void initialize_func_gradient()
            {
                /*set x, v, E and g (this is -grad(E), hence the force) by wrapping position, velocity and gradient arrays in the integrator*/
                _integrator.wrapv(_v); 	//the velocity array wraps the integrator velocity array so that it updates concurrently
                _integrator.wrapf(g_); 	//the force array wraps the integrator force array so that it updates concurrently
                _integrator.wrapx(x_);	//the coordinates array wraps the integrator coordinates array so that it updates concurrently
                _integrator.wrapE(f_);	//the function value (E) wraps the integrator energy so that it updates concurrently
                compute_func_gradient(x_, f_, g_); //compute at initialisation to compute rms_, they'll be recomputed by integrator
                rms_ = norm(g_) / sqrt(N);
                func_initialized_ = true;
            }
  };

  MODIFIED_FIRE::MODIFIEDFIRE(pele::BasePotential * potential, pele::BaseIntegrator * integrator,
		  const pele::Array<double> & x0, double dtstart, double dtmax,
		  size_t Nmin, double finc, double fdec, double fa, double astart, double tol=1e-2):
		  GradientOptimizer(potential,x0,tol=1e-4), //call GradientOptimizer constructor
		  _integrator(integrator(potential, x0, dtstart)),//call to BaseIntegrator constructor to initialise integrator variables
		  _potential(potential), _dtstart(dtstart),
		  _dtmax(dtmax), _astart(astart), _a(astart),
		  _Nmin(Nmin), _finc(finc), _fdec(fdec),
		  _fa(fa), _potential(potential), _dt(dtstart),
    	  _xold(x0.copy()), _N(x0.size())
  	  	  {}

  MODIFIED_FIRE::one_iteration()
  {
	  size_t N = x_.size();
	  nfev_ += 1;
	  iter_number_ += 1;
	  _fire_iter_number += 1; //this is different from iter_number_ which does not get reset

	  _integrator.oneiteration();
	  double P = arraydot(_v,g_);

	  if (P > 0)
	  {
		  /*equation written in this conditional statement _v = (1- _a)*_v + _a * funit * vnorm*/

		  double ifnorm = 1. / norm(g_);
		  double ivnorm = 1. / norm(_v);

		  for (size_t i =0; i < v.size(); ++v)
		  {
			  v[i] += (1. - _a) * _v[i] + _a * g_[i] * _v[i] * ifnorm * ivnorm;
		  }

		  if (_fire_iter_number > _Nmin)
		  {
			  _dt = std::min(_dt* _finc, _dtmax);
			  integrator.set_dt(_dt);
			  _a *= _fa;
		  }
		  _xold(x_.copy());
		  rms_ = norm(g_) / sqrt(N); //update rms
	  }
	  else
	  {
		  _dt *= _fdec;
		  integrator.set_dt(_dt);
		  integrator.set_v_to_vstart(); 		//reset velocity to initial (0)
		  x_(_xold.copy());						//reset position to the one before the step (core of modified fire)
		  _a = _astart;
		  _fire_iter_number = 0;
	  }
  }
}

#endif
