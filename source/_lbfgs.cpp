#include <assert.h>

#include "_lbfgs.h"
#include "base_potential.h"
#include "array.h"
#include <vector>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <iomanip>

using namespace LBFGS_ns;
using std::vector;
using std::cout;

/**
 * compute the dot product of two vectors
 */
double vecdot(std::vector<double> const v1, std::vector<double> const v2)
{
  assert(v1.size() == v2.size());
  size_t i;
  double dot = 0.;
  for (i=0; i<v1.size(); ++i) {
    dot += v1[i] * v2[i];
  }
  return dot;
}

/**
 * compute the L2 norm of a vector
 */
double vecnorm(std::vector<double> const v)
{
  return sqrt(vecdot(v, v));
}

LBFGS::LBFGS(
    pele::BasePotential * potential,
    const pele::Array<double> & x0,
    double tol,
    int M
    )
  :
    potential_(potential),
    M_(M),
    tol_(tol),
    maxstep_(0.1),
    max_f_rise_(1e-4),
    maxiter_(1000),
    iprint_(-1),
    verbosity_(0),
    iter_number_(0),
    nfev_(0),
    H0_(0.1),
    k_(0),
    func_initialized_(false)
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
 * initialize the func and gradient and rms
 */
void LBFGS::initialize_func_gradient()
{
  // compute the func and gradient at the current locations
  // and store them
  size_t N = x_.size();
  compute_func_gradient(x_, f_, g_);
  rms_ = vecnorm(g_) / sqrt(N);
  func_initialized_ = true;
}

/**
 * Set the initial func and gradient.  This can be used
 * to avoid one potential call
 */
void LBFGS::set_func_gradient(double f, pele::Array<double> grad)
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
  }
  rms_ = vecnorm(g_) / sqrt(g_.size());
  func_initialized_ = true;
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

void LBFGS::run(int const niter)
{
    if (! func_initialized_){
      // note: this needs to be both here and in one_iteration
      initialize_func_gradient();
    }

    // iterate until the stop criterion is satisfied or maximum number of
    // iterations is reached
    for (int i = 0; i < niter; ++i)
    {
      if (stop_criterion_satisfied()){
        break;
      }
      one_iteration();
    }

}
void LBFGS::run()
{
  run(maxiter_ - iter_number_);
}

void LBFGS::update_memory(
          std::vector<double> & xold,
          std::vector<double> & gold,
          std::vector<double> & xnew,
          std::vector<double> & gnew
          )
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

bool LBFGS::stop_criterion_satisfied()
{
  return rms_ <= tol_;
}

void LBFGS::compute_func_gradient(std::vector<double> & x, double & func,
      std::vector<double> & gradient)
{
  nfev_ += 1;
  // wrap the vectors as pele::Array objects
  pele::Array<double> xarray(&x[0], x.size());
  pele::Array<double> garray(&gradient[0], gradient.size());

  // pass the arrays to the potential
  func = potential_->get_energy_gradient(x, gradient);
}

void LBFGS::set_H0(double H0)
{
  if (iter_number_ > 0){
    cout << "warning: setting H0 after the first iteration.\n";
  }
  H0_ = H0;
}
