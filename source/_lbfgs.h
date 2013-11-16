#ifndef _PELE_LBFGS_H__
#define _PELE_LBFGS_H__

#include <vector>
#include "base_potential.h"
#include "array.h"
#include "optimizer.h"

using std::vector;

namespace pele{
  /**
   * An implementation of the LBFGS optimization algorithm in c++.  This
   * Implementation uses a backtracking linesearch.
   */
  class LBFGS : public GradientOptimizer{
    private : 

      int M_; /**< The length of the LBFGS memory */
      double max_f_rise_; /**< The maximum the function is allowed to rise in a
                           * given step.  This is the criterion for the
                           * backtracking line search.
                           */

      // places to store the lbfgs memory
      std::vector<vector<double> > s_;
      std::vector<vector<double> > y_;
      std::vector<double> rho_;
      double H0_;
      int k_; /**< Counter for how many times the memory has been updated */

      // 
      std::vector<double> step_;

    public :
      /**
       * Constructor
       */
      LBFGS(
          pele::BasePotential * potential,
          const pele::Array<double> & x0,
          double tol=1e-4,
          int M=4
          );

      /**
       * Destructor
       */
      ~LBFGS() {}

      /**
       * Do one iteration iteration of the optimization algorithm
       */
      void one_iteration();

      // functions for setting the parameters
      void set_H0(double H0)
      {
          if (iter_number_ > 0){
              cout << "warning: setting H0 after the first iteration.\n";
          }
          H0_ = H0;
      }
      void set_max_f_rise(double max_f_rise) { max_f_rise_ = max_f_rise; }

      // functions for accessing the results
      double get_H0() { return H0_; }

    private :

      /**
       * Add a step to the LBFGS Memory
       * This updates s_, y_, rho_, H0_, and k_
       */
      void update_memory(
          std::vector<double> & xold,
          std::vector<double> & gold,
          std::vector<double> & xnew,
          std::vector<double> & gnew
          );

      /**
       * Compute the LBFGS step from the memory
       */
      void compute_lbfgs_step();

      /**
       * Take the step and do a backtracking linesearch if necessary.
       * Apply the maximum step size constraint and ensure that the function
       * does not rise more than the allowed amount.
       */
      double backtracking_linesearch();

  };

  LBFGS::LBFGS(
          pele::BasePotential * potential,
          const pele::Array<double> & x0,
          double tol,
          int M)
  :
          GradientOptimizer(potential, x0, tol),
          M_(M),
          max_f_rise_(1e-4),
          H0_(0.1),
          k_(0)
  {
      // set the precision of the printing
      cout << std::setprecision(12);

      size_t N = x0.size();
      // allocate arrays
      x_ = std::vector<double>(N);
      g_ = std::vector<double>(N);

      y_ = std::vector<vector<double> >(M_, vector<double>(N));
      s_ = std::vector<vector<double> >(M_, vector<double>(N));
      rho_ = std::vector<double>(M_);
      step_ = std::vector<double>(N);

      for (size_t j2 = 0; j2 < N; ++j2){
          x_[j2] = x0[j2];
      }
  }

  /**
   * Do one iteration iteration of the optimization algorithm
   */
  void LBFGS::one_iteration()
  {
      if (! func_initialized_)
          initialize_func_gradient();

      std::vector<double> x_old = x_;
      std::vector<double> g_old = g_;

      compute_lbfgs_step();

      double stepsize = backtracking_linesearch();

      update_memory(x_old, g_old, x_, g_);
      if ((iprint_ > 0) && (iter_number_ % iprint_ == 0)){
          cout << "lbgs: " << iter_number_
                  << " f " << f_
                  << " rms " << rms_
                  << " stepsize " << stepsize << "\n";
      }
      iter_number_ += 1;
  }

  void LBFGS::update_memory(
          std::vector<double> & xold,
          std::vector<double> & gold,
          std::vector<double> & xnew,
          std::vector<double> & gnew)
  {
      // update the lbfgs memory
      // This updates s_, y_, rho_, and H0_, and k_
      int klocal = k_ % M_;
      for (size_t j2 = 0; j2 < x_.size(); ++j2){
          y_[klocal][j2] = gnew[j2] - gold[j2];
          s_[klocal][j2] = xnew[j2] - xold[j2];
      }

      double ys = vecdot(y_[klocal], s_[klocal]);
      if (ys == 0.) {
          // should print a warning here
          if (verbosity_ > 0) {
              cout << "warning: resetting YS to 1.\n";
          }
          ys = 1.;

      }

      rho_[klocal] = 1. / ys;

      double yy = vecdot(y_[klocal], y_[klocal]);
      if (yy == 0.) {
          // should print a warning here
          if (verbosity_ > 0) {
              cout << "warning: resetting YY to 1.\n";
          }
          yy = 1.;
      }
      H0_ = ys / yy;
      //  cout << "    setting H0 " << H0_
      //    << " ys " << ys
      //    << " yy " << yy
      //    << " rho[i] " << rho_[klocal]
      //    << "\n";

      // increment k
      k_ += 1;

  }

  void LBFGS::compute_lbfgs_step()
  {
      if (k_ == 0){
          double gnorm = vecnorm(g_);
          if (gnorm > 1.) gnorm = 1. / gnorm;
          for (size_t j2 = 0; j2 < x_.size(); ++j2){
              step_[j2] = - gnorm * H0_ * g_[j2];
          }
          return;
      }

      step_ = g_;

      int jmin = std::max(0, k_ - M_);
      int jmax = k_;
      int i;
      double beta;
      vector<double> alpha(M_);

      // loop backwards through the memory
      for (int j = jmax - 1; j >= jmin; --j){
          i = j % M_;
          //cout << "    i " << i << " j " << j << "\n";
          alpha[i] = rho_[i] * vecdot(s_[i], step_);
          for (size_t j2 = 0; j2 < step_.size(); ++j2){
              step_[j2] -= alpha[i] * y_[i][j2];
          }
      }

      // scale the step size by H0
      for (size_t j2 = 0; j2 < step_.size(); ++j2){
          step_[j2] *= H0_;
      }

      // loop forwards through the memory
      for (int j = jmin; j < jmax; ++j){
          i = j % M_;
          //cout << "    i " << i << " j " << j << "\n";
          beta = rho_[i] * vecdot(y_[i], step_);
          for (size_t j2 = 0; j2 < step_.size(); ++j2){
              step_[j2] += s_[i][j2] * (alpha[i] - beta);
          }
      }

      // invert the step to point downhill
      for (size_t j2 = 0; j2 < x_.size(); ++j2){
          step_[j2] *= -1;
      }

  }

  double LBFGS::backtracking_linesearch()
  {
      vector<double> xnew(x_.size());
      vector<double> gnew(x_.size());
      double fnew;

      // if the step is pointing uphill, invert it
      if (vecdot(step_, g_) > 0.){
          if (verbosity_ > 1) {
              cout << "warning: step direction was uphill.  inverting\n";
          }
          for (size_t j2 = 0; j2 < step_.size(); ++j2){
              step_[j2] *= -1;
          }
      }

      double factor = 1.;
      double stepsize = vecnorm(step_);

      // make sure the step is no larger than maxstep_
      if (factor * stepsize > maxstep_){
          factor = maxstep_ / stepsize;
      }

      int nred;
      int nred_max = 10;
      for (nred = 0; nred < nred_max; ++nred){
          for (size_t j2 = 0; j2 < xnew.size(); ++j2){
              xnew[j2] = x_[j2] + factor * step_[j2];
          }
          compute_func_gradient(xnew, fnew, gnew);

          double df = fnew - f_;
          if (df < max_f_rise_){
              break;
          } else {
              factor /= 10.;
              if (verbosity_ > 2) {
                  cout
                  << "energy increased: " << df
                  << " reducing step size to " << factor * stepsize
                  << " H0 " << H0_ << "\n";
              }
          }
      }

      if (nred >= nred_max){
          // possibly raise an error here
          if (verbosity_ > 0) {
              cout << "warning: the line search backtracked too many times\n";
          }
      }

      x_ = xnew;
      g_ = gnew;
      f_ = fnew;
      rms_ = vecnorm(gnew) / sqrt(gnew.size());
      return stepsize * factor;
  }

}

#endif
