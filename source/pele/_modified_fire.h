#ifndef _PELE_MODIFIED_FIRE_H__
#define _PELE_MODIFIED_FIRE_H__

#include <vector>
#include "base_potential.h"
#include "array.h"
#include "optimizer.h"
#include "base_integrator.h"
#include "velocity_verlet.h"
#include "forward_euler.h"

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
	  double _dtstart, _dt, _dtmax, _maxstep, _Nmin, _finc, _fdec, _fa, _astart, _a, _fold;
	  pele::Array<double> _v, _xold, _gold;
	  size_t _fire_iter_number;
	  pele::VelocityVerlet _integrator; //create VelocityVerlet integrator
	  //pele::ForwardEuler _integrator; //create ForwardEuler integrator

    public :

	    /**
		* Constructor
		*/
	  MODIFIED_FIRE(pele::BasePotential * potential, pele::Array<double>& x0, double dtstart, double dtmax, double maxstep,
    		  size_t Nmin=5, double finc=1.1, double fdec=0.5, double fa=0.99, double astart=0.1, double tol=1e-3);
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
    	  	  	nfev_ += 1; 					//this accounts for the energy evaluation done by the integrator
    	  	  	_integrator.wrap_v(_v); 		//the velocity array wraps the integrator velocity array so that it updates concurrently
                _integrator.wrap_g(g_); 		//the gradient array wraps the integrator gradient array so that it updates concurrently
                _integrator.wrap_E(f_);			//the function value (E) wraps the integrator energy so that it updates concurrently
                _integrator.wrap_gold(_gold); 	//the gradient array wraps the integrator gradient array so that it updates concurrently
                _fold = f_;
                rms_ = norm(g_) / sqrt(g_.size());
                for(int k=0; k<x_.size();++k) //set initial velocities (using forward Euler)
				  {
					  _v[k] = -g_[k]*_dt;
				  }
                func_initialized_ = true;
            }

      //consider removing set func gradient, it doesn't fit in this wrapping framework
      void set_func_gradient(double f, Array<double> grad)
          {
    	  	  throw std::runtime_error("MODIFIED_FIRE::set_func_gradient: this function is not implemented "
    	  			  "because the gradient is already initialised by BaseIntegrator");
          }
  };

  MODIFIED_FIRE::MODIFIED_FIRE(pele::BasePotential * potential, pele::Array<double>& x0, double dtstart, double dtmax, double maxstep,
		  size_t Nmin, double finc, double fdec, double fa, double astart, double tol):
		  GradientOptimizer(potential,x0,tol=1e-6), //call GradientOptimizer constructor
		  _dtstart(dtstart), _dt(dtstart),
		  _dtmax(dtmax), _maxstep(maxstep), _Nmin(Nmin),
		  _finc(finc), _fdec(fdec), _fa(fa),
		  _astart(astart), _a(astart), _fold(f_),
		  _v(x0.size(),0), _xold(x0.copy()),_gold(g_.copy()),
		  _fire_iter_number(0),
		  _integrator(potential_, x_, _dtstart, _maxstep)
  	  	  {}

  void MODIFIED_FIRE::one_iteration()
  {
	  size_t k;
	  size_t N = x_.size();
	  nfev_ += 1;
	  iter_number_ += 1;
	  _fire_iter_number += 1; //this is different from iter_number_ which does not get reset
	  double P = -1 * dot(_v,g_);

	  if (P > 0)
	  {
		  //save old configuration in case next step has P < 0

		  _fold = f_; //set f_ old before integration step
		  for(k=0; k<x_.size();++k) //set next xold to current (just updated) x
			  {
				  _xold[k] = x_[k];
			  }

		  /*equation written in this conditional statement _v = (1- _a)*_v + _a * funit * vnorm*/

		  double ifnorm = 1. / norm(g_);
		  double vnorm = norm(_v);

		  for (size_t i =0; i < _v.size(); ++i)
		  {
			  _v[i] = (1. - _a) * _v[i] - _a * g_[i] * ifnorm * vnorm;
		  }

		  if (_fire_iter_number > _Nmin)
		  {
			  /*_dt = _integrator.get_dt();  //dt might have been varied in integrator to adjust maxstep*/
			  _dt = std::min(_dt* _finc, _dtmax);
			  _a *= _fa;
			  _integrator.set_dt(_dt);
		  }

		  rms_ = 1. / (ifnorm * sqrt(N)); //update rms
	  }
	  else
	  {
		  /*_integrator.get_dt();  //dt might have been varied in integrator to adjust maxstep*/
		  _dt *= _fdec;
		  _a = _astart;
		  _fire_iter_number = 0;

		  for(k=0; k<x_.size();++k) //reset position and gradient to the one before the step (core of modified fire) reset velocity to initial (0)
		  	  {
			  	  x_[k] = _xold[k];
			  	  g_[k] = _gold[k]; //gold is updated in velocity verlet and its wrapped in function_gradient_initializer
			  	  _v[k] = 0;
			  	  f_ = _fold;
		  	  }
		  _integrator.set_dt(_dt);
	  }
	  _integrator.oneiteration();
  }
}

#endif
