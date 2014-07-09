#ifndef _PELE_MODIFIED_FIRE_H__
#define _PELE_MODIFIED_FIRE_H__

#include <vector>
#include "base_potential.h"
#include "array.h"
#include "optimizer.h"

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
    double _dtstart, _dt, _dtmax, _maxstep, _Nmin, _finc, _fdec, _fa, _astart, _a, _fold, _ifnorm, _vnorm;
    pele::Array<double> _v, _dx, _xold, _gold;
    size_t _fire_iter_number, _N;
    bool _stepback;
    inline void _ForwardEuler_integration();
    inline void _VelocityVerlet_integration();
public :

      /**
      * Constructor
      */
    MODIFIED_FIRE(std::shared_ptr<pele::BasePotential> potential, pele::Array<double>& x0,
            double dtstart, double dtmax, double maxstep, size_t Nmin=5, 
            double finc=1.1, double fdec=0.5, double fa=0.99, double astart=0.1,
            double tol=1e-4, bool stepback=true);
    /**
     * Destructorgit undo rebase
     */
    virtual ~MODIFIED_FIRE() {}

    /**
     * Do one iteration iteration of the optimization algorithm
     */

    void one_iteration();

    /**
     * Overload initialize_func_gradient from parent class
     */

    void initialize_func_gradient();

    void set_func_gradient(double f, Array<double> grad);

    inline void reset(Array<double> &x0);

  };

  inline void MODIFIED_FIRE::reset(Array<double> &x0)
  {
        //arrays are already wrapped by the integrator, must not wrap them again, just update their values, dont's use array.copy()
        // or will wrap a copy to the array
      if (!func_initialized_) {
          initialize_func_gradient();
      }
      iter_number_ = 0;
      x_.assign(x0);
      f_ = potential_->get_energy_gradient(x_, g_);
      nfev_ = 1;
      //fire specific
      _fire_iter_number = 0;
      _dt = _dtstart;
      _a = _astart;
      _fold = f_;
      _xold.assign(x_);
      _gold.assign(g_);
      for(size_t k=0; k<x_.size();++k){
          _v[k] = -g_[k]*_dt;
      }
      _ifnorm = 1./norm(g_);
      _vnorm = norm(_v);
      rms_ = 1. / (_ifnorm*sqrt(_N));
  }

  inline void MODIFIED_FIRE::_VelocityVerlet_integration()
  {
      /* the minuses in the following expressions are due to the fact that
       * the gradients rather than the forces appear in the expression
       */
      for(size_t i=0; i<_N; ++i) {
          _v[i] -= 0.5 * _dt * (_gold[i] + g_[i]);         //update velocity assumes all masses 1
          _dx[i] = _dt * (_v[i] - 0.5 * _dt * g_[i]);      //build displacement vector, assumes all masses 1
      }
      _gold.assign(g_);             //save gradient as old g
      double normdx = norm(_dx);

      if(normdx > _maxstep){
          _dx *= (_maxstep / normdx); //resize displacement vector is greater than _maxstep
      }

      x_ += _dx;

      f_ = potential_->get_energy_gradient(x_, g_);    //update gradient
  }

  inline void MODIFIED_FIRE::_ForwardEuler_integration()
  {
      /* the minuses in the following expressions are due to the fact that
       * the gradients rather than the forces appear in the expression
       */
      for(size_t i=0; i<_N; ++i) { //this was after get_energy_gradient, moved for testing
          _v[i] -= _dt * g_[i];     //update velocity, assumes all masses are 1
          _dx[i] = _dt * _v[i];     //build displacement vector
      }

      _gold.assign(g_);             //save gradient as old g
      double normdx = norm(_dx);

      if(normdx > _maxstep){
          _dx *= (_maxstep / normdx); //resize displacement vector is greater than _maxstep
      }

      x_ += _dx;

      f_ = potential_->get_energy_gradient(x_, g_);    //update gradient
  }

  inline void MODIFIED_FIRE::one_iteration()
  {
      nfev_ += 1;
      iter_number_ += 1;
      _fire_iter_number += 1; //this is different from iter_number_ which does not get reset

      //save old configuration in case next step has P < 0
      _fold = f_; //set f_ old before integration step
      _xold.assign(x_); //save x as xold, gold saved in integrator (because velocity verlet needs it, if vv is used)

      /*equation written in this conditional statement _v = (1- _a)*_v + _a * funit * vnorm*/

      for (size_t i=0; i < _N; ++i) {
          _v[i] = (1. - _a) * _v[i] - _a * g_[i] * _ifnorm * _vnorm;
      }

      /*run MD*/
      this->_ForwardEuler_integration();
      //this->_VelocityVerlet_integration();

      double P = -1 * dot(_v,g_);

      if (P > 0) {
          if (_fire_iter_number > _Nmin) {
              _dt = std::min(_dt* _finc, _dtmax);
              _a *= _fa;
          }

          _ifnorm = 1./norm(g_);
          _vnorm = norm(_v);
          rms_ = 1. / (_ifnorm * sqrt(_N)); //update rms
      } else {
          _dt *= _fdec;
          _a = _astart;
          _fire_iter_number = 0;
          _v.assign(0);

          //reset position and gradient to the one before the step (core of modified fire) reset velocity to initial (0)
          if (_stepback == true) {
              f_ = _fold;
              //reset position and gradient to the one before the step (core of modified fire) reset velocity to initial (0)
              x_.assign(_xold);
              g_.assign(_gold);
          }
      }
  }
}

#endif
