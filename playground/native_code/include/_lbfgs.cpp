#include "_lbfgs.h"
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
double vecdot(std::vector<double> const v1, std::vector<double> const v2,
    size_t N)
{
  size_t i;
  double dot = 0.;
  for (i=0; i<N; ++i) {
    dot += v1[i] * v2[i];
  }
  return dot;
}

/**
 * compute the L2 norm of a vector
 */
double vecnorm(std::vector<double> const v, size_t N)
{
  return sqrt(vecdot(v, v, N));
}

LBFGS::LBFGS(
    double (*func)(double *, double *, int), 
    double const * x0, 
    int N, 
    int M
    //double tol,
    //double maxstep,
    //double max_f_rise,
    //double H0,
    //int maxiter
    )
  :
    func_f_grad_(func),
    N_(N),
    M_(M),
    tol_(1e-4),
    maxstep_(0.2),
    max_f_rise_(1e-4),
    maxiter_(1000),
    iprint_(-1),
    iter_number_(0),
    nfev_(0),
    H0_(0.1),
    k_(0)
{
  // set the precision of the printing
  cout << std::setprecision(12);

  // allocate arrays
  x_ = std::vector<double>(N_);
  g_ = std::vector<double>(N_);

  y_ = std::vector<vector<double> >(M_, vector<double>(N_));
  s_ = std::vector<vector<double> >(M_, vector<double>(N_));
  rho_ = std::vector<double>(M_);
  step_ = std::vector<double>(N_);

  // set up the current location
  for (int j2 = 0; j2 < N_; ++j2){ 
    x_[j2] = x0[j2];
  }
  compute_func_gradient(x_, f_, g_);
  rms_ = vecnorm(g_, N_) / sqrt(N_);
}

LBFGS::~LBFGS()
{
}

void LBFGS::one_iteration()
{
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

void LBFGS::run()
{
  // iterate until the stop criterion is satisfied or maximum number of
  // iterations is reached
  while (iter_number_ < maxiter_)
  {
    if (stop_criterion_satisfied()){
      break;
    }
    one_iteration();
  }
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
  for (int j2 = 0; j2 < N_; ++j2){ 
    y_[klocal][j2] = gnew[j2] - gold[j2];
    s_[klocal][j2] = xnew[j2] - xold[j2];
  }

  double ys = vecdot(y_[klocal], s_[klocal], N_);
  if (ys == 0.) {
    // should print a warning here
    cout << "warning: resetting YS to 1.\n";
    ys = 1.;
  }

  rho_[klocal] = 1. / ys;

  double yy = vecdot(y_[klocal], y_[klocal], N_);
  if (yy == 0.) {
    // should print a warning here
    cout << "warning: resetting YY to 1.\n";
    yy = 1.;
  }
  H0_ = ys / yy;
//  cout << "    setting H0 " << H0_ 
//    << " ys " << ys 
//    << " yy " << yy 
//    << " rho[i] " << rho_[klocal] 
//    << "\n";

  // incriment k
  k_ += 1;
  
}

void LBFGS::compute_lbfgs_step()
{
  if (k_ == 0){ 
    double gnorm = vecnorm(g_, N_);
    if (gnorm > 1.) gnorm = 1. / gnorm;
    for (int j2 = 0; j2 < N_; ++j2){
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
    alpha[i] = rho_[i] * vecdot(s_[i], step_, N_);
    for (int j2 = 0; j2 < N_; ++j2){
      step_[j2] -= alpha[i] * y_[i][j2];
    }
  }

  // scale the step size by H0
  for (int j2 = 0; j2 < N_; ++j2){
    step_[j2] *= H0_;
  }

  // loop forwards through the memory
  for (int j = jmin; j < jmax; ++j){
    i = j % M_;
    //cout << "    i " << i << " j " << j << "\n";
    beta = rho_[i] * vecdot(y_[i], step_, N_);
    for (int j2 = 0; j2 < N_; ++j2){
      step_[j2] += s_[i][j2] * (alpha[i] - beta);
    }
  }

  // invert the step to point downhill
  for (int j2 = 0; j2 < N_; ++j2){
    step_[j2] *= -1;
  }

}

double LBFGS::backtracking_linesearch()
{
  vector<double> xnew(N_);
  vector<double> gnew(N_);
  double fnew;

  // if the step is pointing uphill, invert it
  if (vecdot(step_, g_, N_) > 0.){
    cout << "warning: step direction was uphill.  inverting\n";
    for (int j2 = 0; j2 < N_; ++j2){
      step_[j2] *= -1;
    }
  }

  double factor = 1.;
  double stepsize = vecnorm(step_, N_);

  // make sure the step is no larger than maxstep_
  if (factor * stepsize > maxstep_){
    factor = maxstep_ / stepsize;
  }

  int nred;
  int nred_max = 10;
  for (nred = 0; nred < nred_max; ++nred){
    for (int j2 = 0; j2 < N_; ++j2){
      xnew[j2] = x_[j2] + factor * step_[j2];
    }
    compute_func_gradient(xnew, fnew, gnew);

    double df = fnew - f_;
    if (df < max_f_rise_){
      break;
    } else {
      factor /= 10.;
      cout 
        << "energy increased: " << df 
        << " reducing step size to " << factor * stepsize 
        << " H0 " << H0_ << "\n";
    }
  }

  if (nred >= nred_max){
    // possibly raise an error here
    cout << "warning: the line search backtracked too many times\n";
  }

  x_ = xnew;
  g_ = gnew;
  f_ = fnew;
  rms_ = vecnorm(gnew, N_) / sqrt(N_);
  return stepsize * factor;
}

int LBFGS::stop_criterion_satisfied()
{
  return rms_ <= tol_;
}

void LBFGS::compute_func_gradient(std::vector<double> & x, double & func,
      std::vector<double> & gradient)
{
  nfev_ += 1;
  func = (*func_f_grad_)(&x[0], &gradient[0], N_);
}

void LBFGS::set_H0(double H0)
{
  if (iter_number_ > 0){
    cout << "warning: setting H0 after the first iteration.\n";
  }
  H0_ = H0;
}
