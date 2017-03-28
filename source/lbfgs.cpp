#include "pele/lbfgs.h"
#include <memory>

using std::cout;

namespace pele {

LBFGS::LBFGS( std::shared_ptr<pele::BasePotential> potential, const pele::Array<double> x0,
        double tol, int M)
    : GradientOptimizer(potential, x0, tol),
      M_(M),
      max_f_rise_(1e-4),
      use_relative_f_(false),
      rho_(M_),
      H0_(0.1),
      k_(0),
      alpha(M_),
      xold(x_.size()),
      gold(x_.size()),
      step(x_.size())
{
    // set the precision of the printing
    cout << std::setprecision(12);

    inv_sqrt_size = 1 / sqrt(x_.size());

    // allocate space for s_ and y_
    s_ = Array<double>(x_.size() * M_);
    y_ = Array<double>(x_.size() * M_);
}

/**
* Do one iteration iteration of the optimization algorithm
*/
void LBFGS::one_iteration()
{
    if (!func_initialized_) {
        initialize_func_gradient();
    }

    // make a copy of the position and gradient
    xold.assign(x_);
    gold.assign(g_);

    // get the stepsize and direction from the LBFGS algorithm
    compute_lbfgs_step(step);

    // reduce the stepsize if necessary
    double stepsize = backtracking_linesearch(step);

    // update the LBFGS memeory
    update_memory(xold, gold, x_, g_);

    // print some status information
    if ((iprint_ > 0) && (iter_number_ % iprint_ == 0)){
        cout << "lbgs: " << iter_number_
            << " E " << f_
            << " rms " << rms_
            << " nfev " << nfev_
            << " stepsize " << stepsize << "\n";
    }
    iter_number_ += 1;
}

void LBFGS::update_memory(
        Array<double> x_old,
        Array<double> g_old,
        Array<double> x_new,
        Array<double> g_new)
{
    // update the lbfgs memory
    // This updates s_, y_, rho_, and H0_, and k_
    int klocal = k_ % M_;
    double ys = 0;
    double yy = 0;
    #pragma simd reduction( + : ys, yy)
    for (size_t j2 = 0; j2 < x_.size(); ++j2){
        size_t ind_j2 = klocal * x_.size() + j2;
        y_[ind_j2] = g_new[j2] - g_old[j2];
        s_[ind_j2] = x_new[j2] - x_old[j2];
        ys += y_[ind_j2] * s_[ind_j2];
        yy += y_[ind_j2] * y_[ind_j2];
    }

    if (ys == 0.) {
        if (verbosity_ > 0) {
            cout << "warning: resetting YS to 1.\n";
        }
        ys = 1.;
    }

    rho_[klocal] = 1. / ys;

    if (yy == 0.) {
        if (verbosity_ > 0) {
            cout << "warning: resetting YY to 1.\n";
        }
        yy = 1.;
    }
    H0_ = ys / yy;

    // increment k
    k_ += 1;
}

void LBFGS::compute_lbfgs_step(Array<double> step)
{
    if (k_ == 0){
        // take a conservative first step
        double gnorm = norm(g_);
        if (gnorm > 1.) {
            gnorm = 1. / gnorm;
        }
        double prefactor =  -gnorm * H0_;
        #pragma simd
        for (size_t j2 = 0; j2 < x_.size(); ++j2){
            step[j2] = prefactor * g_[j2];
        }
        return;
    }

    // copy the gradient into step
    step.assign(g_);

    int jmin = std::max(0, k_ - M_);
    int jmax = k_;
    int i;

    alpha.assign(0.0);
    // loop backwards through the memory
    for (int j = jmax - 1; j >= jmin; --j) {
        i = j % M_;
        double alpha_tmp = 0;
        #pragma simd reduction(+ : alpha_tmp)
        for (size_t j2 = 0; j2 < step.size(); ++j2){
            alpha_tmp += rho_[i] * s_[i * step.size() + j2] * step[j2];
        }
        #pragma simd
        for (size_t j2 = 0; j2 < step.size(); ++j2){
            step[j2] -= alpha_tmp * y_[i * step.size() + j2];
        }
        alpha[i] = alpha_tmp;
    }

    // scale the step size by H0, invert the step to point downhill
    #pragma simd
    for (size_t j2 = 0; j2 < step.size(); ++j2){
        step[j2] *= -H0_;
    }

    // loop forwards through the memory
    for (int j = jmin; j < jmax; ++j) {
        i = j % M_;
        double beta = 0;
        #pragma simd reduction(+ : beta)
        for (size_t j2 = 0; j2 < step.size(); ++j2){
            beta -= rho_[i] * y_[i * step.size() + j2] * step[j2];  // -= due to inverted step
        }
        double alpha_beta = alpha[i] - beta;
        #pragma simd
        for (size_t j2 = 0; j2 < step.size(); ++j2){
            step[j2] -= s_[i * step.size() + j2] * alpha_beta;  // -= due to inverted step
        }
    }
}

double LBFGS::backtracking_linesearch(Array<double> step)
{
    double fnew;

    // if the step is pointing uphill, invert it
    if (dot(step, g_) > 0.){
        if (verbosity_ > 1) {
            cout << "warning: step direction was uphill.  inverting\n";
        }
        #pragma simd
        for (size_t j2 = 0; j2 < step.size(); ++j2){
            step[j2] *= -1;
        }
    }

    double factor = 1.;
    double stepsize = norm(step);

    // make sure the step is no larger than maxstep_
    if (factor * stepsize > maxstep_){
        factor = maxstep_ / stepsize;
    }

    int nred;
    int nred_max = 10;
    for (nred = 0; nred < nred_max; ++nred){
        #pragma simd
        for (size_t j2 = 0; j2 < x_.size(); ++j2){
            x_[j2] = xold[j2] + factor * step[j2];
        }
        compute_func_gradient(x_, fnew, g_);

        double df = fnew - f_;
        if (use_relative_f_) {
            double absf = 1e-100;
            if (f_ != 0) {
                absf = std::abs(f_);
            }
            df /= absf;
        }
        if (df < max_f_rise_){
            break;
        }
        else {
            factor *= 0.1;
            if (verbosity_ > 2) {
                cout
                    << "energy increased by " << df
                    << " to " << fnew
                    << " from " << f_
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

    f_ = fnew;
    rms_ = norm(g_) * inv_sqrt_size;
    return stepsize * factor;
}

void LBFGS::reset(pele::Array<double> &x0)
{
    if (x0.size() != x_.size()){
        throw std::invalid_argument("The number of degrees of freedom (x0.size()) cannot change when calling reset()");
    }
    k_ = 0;
    iter_number_ = 0;
    nfev_ = 0;
    x_.assign(x0);
    initialize_func_gradient();
}


}
