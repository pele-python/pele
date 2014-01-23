#ifndef _PELE_MODIFIED_FIRE_H__
#define _PELE_MODIFIED_FIRE_H__

#include <vector>
#include "base_potential.h"
#include "array.h"
#include "optimizer.h"
#include "base_integrator.h"
#include "velocity_verlet.h"

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
	  double _dtstart, _dt, _dtmax, _Nmin, _finc, _fdec, _fa, _astart, _a;
	  pele::Array<double> _v, _xold, _gold;
	  size_t _fire_iter_number, _N;
	  pele::VelocityVerlet _integrator; //create VelocityVerlet integrator

    public :

	    /**
		* Constructor
		*/
	  MODIFIED_FIRE(pele::BasePotential * potential, pele::Array<double> x0, double dtstart, double dtmax,
    		  size_t Nmin=5, double finc=1.1, double fdec=0.5, double fa=0.99, double astart=0.1, double tol=1e-2);
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
    	  	  	_fire_iter_number = 0;
    	  	  	_N = 0;
    	  	  	/*set x, v, E and g (this is grad(E), hence -force) by wrapping position, velocity and gradient arrays in the integrator*/
    	  	  	//x is alredy wrapped by the integrator constructor

    	  	  	compute_func_gradient(x_, f_, g_); //compute at initialisation to compute rms_, they'll be recomputed by integrator
    	  	  	_integrator.wrap_v(_v); 		//the velocity array is wrapped by the integrator velocity array so that it updates concurrently
                _integrator.wrap_g(g_); 		//the gradient array is wrapped by the integrator gradient array so that it updates concurrently
                _integrator.wrap_E(f_);			//the function value (E) wraps the integrator energy so that it updates concurrently
                _integrator.wrap_gold(_gold); 	//the gradient array is wrapped by  the integrator gradient array so that it updates concurrently

                rms_ = norm(g_) / sqrt(g_.size());
                func_initialized_ = true;
            }

      void set_func_gradient(double f, Array<double> grad)
          {
              if (grad.size() != g_.size()){
                  throw std::invalid_argument("the gradient has the wrong size");
              }
              if (iter_number_ > 0){
                  cout << "warning: setting f and grad after the first iteration.  this is dangerous.\n";
              }
              // copy the function and gradient
              f_ = f;
              size_t N = x_.size();
              for (size_t j2 = 0; j2 < N; ++j2){
                  g_[j2] = grad[j2];
                  _gold[j2] = grad[j2];
              }

              rms_ = norm(g_) / sqrt(g_.size());

              //fire specific directives
              _fire_iter_number = 0;
			  _N = 0;

			  _integrator.wrap_v(_v); 		//the velocity array wraps the integrator velocity array so that it updates concurrently
			  _integrator.wrap_g(g_); 		//the gradient array wraps the integrator gradient array so that it updates concurrently
			  _integrator.wrap_E(f_);		//the function value (E) wraps the integrator energy so that it updates concurrently
			  _integrator.wrap_gold(_gold); //the gradient array wraps the integrator gradient array so that it updates concurrently

			  func_initialized_ = true;
          }
  };

  MODIFIED_FIRE::MODIFIED_FIRE(pele::BasePotential * potential, pele::Array<double> x0, double dtstart, double dtmax,
		  size_t Nmin, double finc, double fdec, double fa, double astart, double tol):
		  GradientOptimizer(potential,x0,tol=1e-2), //call GradientOptimizer constructor
		  _dtstart(dtstart), _dt(dtstart),
		  _dtmax(dtmax), _Nmin(Nmin),
		  _finc(finc), _fdec(fdec), _fa(fa),
		  _astart(astart), _a(astart), _v(x0.size(),0),
		  _xold(x0.copy()),_gold(g_.copy()),
		  _N(x0.size()),
		  _integrator(potential_, x_, _dtstart)
  	  	  {}

  void MODIFIED_FIRE::one_iteration()
  {
	  size_t k;
	  size_t N = x_.size();
	  nfev_ += 1;
	  iter_number_ += 1;
	  _fire_iter_number += 1; //this is different from iter_number_ which does not get reset

	  _integrator.oneiteration();
	  double P = -1 * dot(_v,g_);

	  if (P > 0)
	  {
		  /*equation written in this conditional statement _v = (1- _a)*_v + _a * funit * vnorm*/

		  double ifnorm = 1. / norm(g_);
		  double ivnorm = 1. / norm(_v);

		  for (size_t i =0; i < _v.size(); ++i)
		  {
			  _v[i] = (1. - _a) * _v[i] - _a * g_[i] * _v[i] * ifnorm * ivnorm;
		  }

		  if (_fire_iter_number > _Nmin)
		  {
			  _dt = std::min(_dt* _finc, _dtmax);
			  _a *= _fa;
		  }

		  for(k=0; k<x_.size();++k) //set next xold to current x
		  	  {
			  	  _xold[k] = x_[k];
		  	  }
		  rms_ = norm(g_) / sqrt(N); //update rms
	  }
	  else
	  {
		  _dt *= _fdec;
		  _a = _astart;
		  _fire_iter_number = 0;

		  for(k=0; k<x_.size();++k) //reset position and gradient to the one before the step (core of modified fire) reset velocity to initial (0)
		  	  {
			  	  x_[k] = _xold[k];
			  	  g_[k] = _gold[k]; //gold is updated in velocity verlet and its wrapped in function_gradient_initializer
			  	  _v[k] = 0;
		  	  }
	  }

	  _integrator.set_dt(_dt);
  }
}

#endif
