#ifndef _PELE_MODIFIED_FIRE_H__
#define _PELE_MODIFIED_FIRE_H__

#include <vector>
#include "base_potential.h"
#include "array.h"
#include "optimizer.h"
#include "base_integrator.h"
#include "forward_euler_fire.h"

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
   *
   * This implementation of the algorithm differs significantly from the original
   * algorithm in the order in which the steps are taken. Here we do the following:
   *
   * -> initialise the velocities and gradients for some given coordinates
   * -> set the velocity to the fire *inertial* velocity
   *
   * Only then we use the MD integrator that here does things in this order:
   *
   * -> compute the velocity difference and add it to the current velocity
   * -> compute the new position given this velocity
   * -> recompute gradient and energy
   *
   * Once the integrator is done we continue following FIRE and compute
   *
   * P = -g * f
   *
   * here comes the other modification to the algorithm: if stepback is set to false
   * then proceed just as the original FIRE algorithm prescribes, if stepback is set
   * to true (default) then whenever P<=0 we undo the last step besides carrying out
   * the operations defined by the original FIRE algorithm.
   *
   */

    class MODIFIED_FIRE : public GradientOptimizer{
    private :
      double _dtstart, _dt, _dtmax, _maxstep, _Nmin, _finc, _fdec, _fa, _astart, _a, _fold;
      pele::Array<double> _v, _xold, _gold;
      size_t _fire_iter_number, _N;
      bool _stepback;
      //pele::VelocityVerlet _integrator; //create VelocityVerlet integrator
      pele::ForwardEuler _integrator; //create ForwardEuler integrator

    public :

        /**
        * Constructor
        */
      MODIFIED_FIRE(pele::BasePotential * potential, pele::Array<double>& x0, double dtstart, double dtmax, double maxstep,
              size_t Nmin=5, double finc=1.1, double fdec=0.5, double fa=0.99, double astart=0.1, double tol=1e-4, bool stepback=true);
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
                    nfev_ += 1;                     //this accounts for the energy evaluation done by the integrator (g is computed by the constructor of the integrator)
                    _integrator.wrap_v(_v);         //the velocity array wraps the integrator velocity array so that it updates concurrently
                _integrator.wrap_g(g_);         //the gradient array wraps the integrator gradient array so that it updates concurrently
                _integrator.wrap_E(f_);            //the function value (E) wraps the integrator energy so that it updates concurrently
                _integrator.wrap_gold(_gold);     //the gradient array wraps the integrator gradient array so that it updates concurrently
                _fold = f_;
                rms_ = norm(g_) / sqrt(g_.size());
                for(size_t k=0; k<x_.size();++k) //set initial velocities (using forward Euler)
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

      inline void reset(Array<double> &x0)
      {
          //arrays are already wrapped by the integrator, must not wrap them again, just update their values, dont's use array.copy()
          // or will wrap a copy to the array
          iter_number_ = 0;
          x_.assign(x0);
          f_ = potential_->get_energy_gradient(x_, g_);
          nfev_ = 1;
          //fire specific
          _fire_iter_number = 0;
          _dt = _dtstart;
          _a = _astart;
          _fold = f_;
          rms_ = norm(g_) / sqrt(g_.size());
          _xold.assign(x_);
          _gold.assign(g_);
          for(size_t k=0; k<x_.size();++k){
              _v[k] = -g_[k]*_dt;
          }
          _integrator.set_dt(_dt);
      }


    };

  MODIFIED_FIRE::MODIFIED_FIRE(pele::BasePotential * potential, pele::Array<double>& x0, double dtstart, double dtmax, double maxstep,
          size_t Nmin, double finc, double fdec, double fa, double astart, double tol, bool stepback):
          GradientOptimizer(potential,x0,tol=1e-4), //call GradientOptimizer constructor
          _dtstart(dtstart), _dt(dtstart),
          _dtmax(dtmax), _maxstep(maxstep), _Nmin(Nmin),
          _finc(finc), _fdec(fdec), _fa(fa),
          _astart(astart), _a(astart), _fold(f_),
          _v(x0.size(),0), _xold(x0.copy()),_gold(g_.copy()),
          _fire_iter_number(0), _N(x_.size()),
          _stepback(stepback),
          _integrator(potential_, x_, _dtstart, _maxstep)
              {}

  inline void MODIFIED_FIRE::one_iteration()
  {
      double ifnorm, vnorm, P;
      //size_t k;
      nfev_ += 1;
      iter_number_ += 1;
      _fire_iter_number += 1; //this is different from iter_number_ which does not get reset

      //save old configuration in case next step has P < 0
      _fold = f_; //set f_ old before integration step
      _xold.assign(x_); //save x as xold (gold is saved in the integrator)

      /*equation written in this conditional statement _v = (1- _a)*_v + _a * funit * vnorm*/
      ifnorm = 1. / norm(g_);
      vnorm = norm(_v);

      for (size_t i =0; i < _v.size(); ++i)
      {
          _v[i] = (1. - _a) * _v[i] - _a * g_[i] * ifnorm * vnorm;
      }

      /*run MD*/
      _integrator.oneiteration();

      P = -1 * dot(_v,g_);

      if (P > 0)
      {
          if (_fire_iter_number > _Nmin)
          {
              _dt = std::min(_dt* _finc, _dtmax);
              _a *= _fa;
              _integrator.set_dt(_dt);
          }

          rms_ = 1. / (ifnorm * sqrt(_N)); //update rms
      }
      else
      {
          _dt *= _fdec;
          _a = _astart;
          _fire_iter_number = 0;
          _integrator.set_dt(_dt);
          _v.assign(0);

          //reset position and gradient to the one before the step (core of modified fire) reset velocity to initial (0)
          if (_stepback == true)
          {
              f_ = _fold;
              //reset position and gradient to the one before the step (core of modified fire) reset velocity to initial (0)
              x_.assign(_xold);
              g_.assign(_gold);
          }
      }
  }
}

#endif
