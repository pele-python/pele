#ifndef _PELE_OPTIMIZER_H__
#define _PELE_OPTIMIZER_H__

#include "base_potential.h"
#include "array.h"
#include <vector>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <iomanip>

using std::vector;
using std::runtime_error;
using std::cout;

namespace pele{

/**
 * this defines the basic interface for optimizers.  All pele optimizers
 * should derive from this class.
 */
class Optimizer{
public:
    virtual void one_iteration() = 0;

    /**
     * Run the optimization algorithm until the stop criterion is satisfied or
     * until the maximum number of iterations is reached
     */
    virtual void run() = 0;

    /**
     * Run the optimization algorithm for niter iterations or until the
     * stop criterion is satisfied
     */
    virtual void run(int const niter) = 0;

    /**
     * accessors
     */
    virtual Array<double> get_x() = 0;
    virtual Array<double> get_g() = 0;
    virtual double get_f() = 0;
    virtual double get_rms() = 0;
    virtual int get_nfev() = 0;
    virtual int get_niter() = 0;
    virtual bool success() = 0;
};

/**
 * This defines the basic interface for optimizers.  All pele optimizers
 * should derive from this class.
 */
class GradientOptimizer : public Optimizer{
protected :
    // input parameters
    /**
     * A pointer to the object that computes the function and gradient
     */
    pele::BasePotential * potential_;

    double tol_; /**< The tolerance for the rms gradient */
    double maxstep_; /**< The maximum step size */

    int maxiter_; /**< The maximum number of iterations */
    int iprint_; /**< how often to print status information */
    int verbosity_; /**< How much information to print */

    int iter_number_; /**< The current iteration number */
    int nfev_; /**< The number of function evaluations */

    // variables representing the state of the system
    Array<double> x_; /**< The current coordinates */
    double f_; /**< The current function value */
    Array<double> g_; /**< The current gradient */
    double rms_; /**< The root mean square of the gradient */

    /**
     * This flag keeps track of whether the function and gradient have been
     * initialized.  This allows the initial function and gradient to be computed
     * outside of the constructor and also allows the function and gradient to
     * be passed rather than computed.  The downside is that it complicates the
     * logic because this flag must be checked at all places where the gradient,
     * function value, or rms can be first accessed.
     */
    bool func_initialized_;

public :
    GradientOptimizer(pele::BasePotential * potential,
          const pele::Array<double> x0,
          double tol=1e-4)
    :
      potential_(potential),
      tol_(tol),
      maxstep_(0.1),
      maxiter_(1000),
      iprint_(-1),
      verbosity_(0),
      iter_number_(0),
      nfev_(0),
      x_(x0.copy()),
      f_(0.),
      g_(x0.size()),
      rms_(1e10),
      func_initialized_(false)
    {}

    /**
     * Do one iteration iteration of the optimization algorithm
     */
    virtual void one_iteration() = 0;

    /**
     * Run the optimization algorithm until the stop criterion is satisfied or
     * until the maximum number of iterations is reached
     */
    void run(int const niter)
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

    /**
     * Run the optimzation algorithm for niter iterations or until the
     * stop criterion is satisfied
     */
    void run()
    {
        run(maxiter_ - iter_number_);
    }

    /**
     * Set the initial func and gradient.  This can be used
     * to avoid one potential call
     */
    virtual void set_func_gradient(double f, Array<double> grad)
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
        rms_ = norm(g_) / sqrt(g_.size());
        func_initialized_ = true;
    }

    // functions for setting the parameters
    void set_tol(double tol) { tol_ = tol; }
    void set_maxstep(double maxstep) { maxstep_ = maxstep; }
    void set_max_iter(int max_iter) { maxiter_ = max_iter; }
    void set_iprint(int iprint) { iprint_ = iprint; }
    void set_verbosity(int verbosity) { verbosity_ = verbosity; }


    // functions for accessing the status of the optimizer
    Array<double> get_x() { return x_; }
    Array<double> get_g() { return g_; }
    double get_f() { return f_; }
    double get_rms() { return rms_; }
    int get_nfev() { return nfev_; }
    int get_niter() { return iter_number_; }
    int get_maxiter() { return maxiter_; }
    bool success() { return stop_criterion_satisfied(); }

    /**
     * Return true if the termination condition is satisfied, false otherwise
     */
    bool stop_criterion_satisfied()
    {
        if (! func_initialized_) initialize_func_gradient();
        return rms_ <= tol_;
    }

protected :

    /**
     * Compute the func and gradient of the objective function
     */
    void compute_func_gradient(Array<double> x, double & func,
            Array<double> gradient)
    {
        nfev_ += 1;

        // pass the arrays to the potential
        func = potential_->get_energy_gradient(x, gradient);
    }

    /**
     * compute the initial func and gradient
     */
    virtual void initialize_func_gradient()
    {
        // compute the func and gradient at the current locations
        // and store them
        size_t N = x_.size();
        compute_func_gradient(x_, f_, g_);
        rms_ = norm(g_) / sqrt(N);
        func_initialized_ = true;
    }

};
}

#endif
