#ifndef _PELE_OPTIMIZER_H__
#define _PELE_OPTIMIZER_H__

#include <Eigen/Dense>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <vector>

#include "array.hpp"
#include "base_potential.hpp"
#include "preprocessor_directives.hpp"
#include "xsum.h"

namespace pele {

enum StopCriterionType {
  GRADIENT = 0,
  STEPNORM = 1,
  NEWTON = 2,
};

class GradientOptimizer;

class AbstractStopCriterion {
 public:
  virtual bool stop_criterion_satisfied() = 0;
  virtual ~AbstractStopCriterion() = default;
};

class GradientStopCriterion : public AbstractStopCriterion {
 private:
  double tol_;
  GradientOptimizer *opt_;

 public:
  GradientStopCriterion(double tol, GradientOptimizer *opt)
      : tol_(tol), opt_(opt) {}

  bool stop_criterion_satisfied();
};

/**
 * @brief Stop criterion based on the step norm.
 *
 *
 */
class StepnormStopCriterion : public AbstractStopCriterion {
 private:
  double gradient_tol_;
  double step_norm_tol_;
  GradientOptimizer *opt_;

 public:
  StepnormStopCriterion(double tol, GradientOptimizer *opt)
      : gradient_tol_(tol * 1e-3), step_norm_tol_(tol), opt_(opt) {}
  bool stop_criterion_satisfied();
};

class NewtonStopCriterion : public AbstractStopCriterion {
 private:
  double gradient_tol_;
  double newton_tol_;
  double offset_factor_;
  GradientOptimizer *opt_;
  Eigen::MatrixXd hessian_;
  Eigen::VectorXd gradient;
  Eigen::VectorXd eigenvalues_;
  Eigen::VectorXd newton_step_;
  Eigen::MatrixXd eigenvectors_;
  std::shared_ptr<BasePotential> potential_;

 public:
  NewtonStopCriterion(double tol, GradientOptimizer *opt);
  bool stop_criterion_satisfied();
};

// class blah
// {
// public:
//     blah(double a) {std::cout << "this works" << "\n";};
//     virtual ~blah() {};
// };

/**
 * this defines the basic interface for optimizers.  All pele optimizers
 * should derive from this class.
 */
class Optimizer {
 public:
  /**
   * virtual destructor
   */
  virtual ~Optimizer() {}

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
  inline virtual Array<double> get_x() const = 0;
  inline virtual Array<double> get_g() const = 0;
  inline virtual double get_f() const = 0;
  inline virtual double get_rms() const = 0;
  inline virtual int get_nfev() const = 0;
  inline virtual int get_niter() const = 0;
  inline virtual bool success() = 0;
  inline virtual int get_nhev() const = 0;
};

/**
 * This defines the basic interface for optimizers.  All pele optimizers
 * should derive from this class.
 */
class GradientOptimizer : public Optimizer {
 protected:
  // input parameters
  /**
   * A pointer to the object that computes the function and gradient
   */
  std::shared_ptr<pele::BasePotential> potential_;
  // friend class LineSearch;  // LineSearch can access internals
  double tol_;     /**< The tolerance for the rms gradient */
  double maxstep_; /**< The maximum step size */
  int maxiter_;    /**< The maximum number of iterations */
  int iprint_;     /**< how often to print status information */
  int verbosity_;  /**< How much information to print */

  int iter_number_; /**< The current iteration number */
  /**
   * number of get_energy_gradient evaluations.
   */
  int nfev_;
  /**
   * number of get_energy_gradient_hessian evaluations
   */
  int nhev_;

  // variables representing the state of the system
  Array<double> x_;      /**< The current coordinates */
  double f_;             /**< The current function value */
  Array<double> g_;      /**< The current gradient */
  double gradient_norm_; /**< The root mean square of the gradient */
  double step_norm_;     /**< The norm of the step taken */
  std::unique_ptr<AbstractStopCriterion> stop_criterion;

  /**
   * This flag keeps track of whether the function and gradient have been
   * initialized.This allows the initial function and gradient to be
   * computed outside of the constructor and also allows the function and
   * gradient to be passed rather than computed.  The downside is that it
   * complicates the logic because this flag must be checked at all places
   * where the gradient, function value, or rms can be first accessed.
   */
  bool func_initialized_;

  bool last_step_failed_; /**< Whether the last step failed
  helpful for newton, when combining in mixed descent */

  bool succeeded_; /**< Whether the optimization succeeded */

  bool save_trajectory_; /**< Whether to save the trajectory */

  int iterations_before_save_; /**< How often to save the trajectory */

 public:
  GradientOptimizer(std::shared_ptr<pele::BasePotential> potential,
                    const pele::Array<double> x0, double tol = 1e-4,
                    bool save_trajectory = false,
                    int iterations_before_save = 1,
                    StopCriterionType stop = GRADIENT)
      : potential_(potential),
        tol_(tol),
        maxstep_(0.1),
        maxiter_(1000),
        iprint_(-1),
        verbosity_(0),
        iter_number_(0),
        nfev_(0),
        nhev_(0),
        x_(x0.copy()),
        f_(0.),
        g_(x0.size()),
        gradient_norm_(1e10),
        step_norm_(0),
        func_initialized_(false),
        last_step_failed_(false),
        succeeded_(false),
        save_trajectory_(save_trajectory),
        iterations_before_save_(iterations_before_save) {
    if (stop == GRADIENT) {
      stop_criterion = std::make_unique<GradientStopCriterion>(tol_, this);
    } else if (stop == STEPNORM) {
      stop_criterion = std::make_unique<StepnormStopCriterion>(tol_, this);
    } else if (stop == NEWTON) {
      stop_criterion = std::make_unique<NewtonStopCriterion>(tol_, this);
    } else {
      throw std::invalid_argument("invalid stop criterion");
    }
  }

  virtual ~GradientOptimizer() {}

  // copy constructor raises exception that it's not implemented
  GradientOptimizer(const GradientOptimizer &) : Optimizer() {
    throw std::runtime_error(
        "GradientOptimizer copy constructor not implemented");
  }

  GradientOptimizer() : Optimizer() {
    throw std::runtime_error(
        "GradientOptimizer default constructor not "
        "implemented. Check derived class function calls");
  }

  /**
   * Do one iteration iteration of the optimization algorithm
   */
  virtual void one_iteration() = 0;

  /**
   * Run the optimization algorithm until the stop criterion is satisfied or
   * until the maximum number of iterations is reached
   */
  virtual void run(int const niter) {
    if (!func_initialized_) {
      // note: this needs to be both here and in one_iteration
      initialize_func_gradient();
    }

    // iterate until the stop criterion is satisfied or maximum number of
    // iterations is reached
    for (int i = 0; i < niter; ++i) {
      if (stop_criterion_satisfied()) break;
      one_iteration();
      if (save_trajectory_) {
        // costly info takes up more memory, so only save when necessary
        if (iter_number_ % iterations_before_save_ == 0)
          update_costly_trajectory_info();
        // cheap info takes up less memory, so save every iteration
        update_cheap_trajectory_info();
      }
    }
  }

  void print_stats() {
    std::cout << "iter: " << iter_number_ << std::endl;
    std::cout << " rms: " << gradient_norm_ << std::endl;
    std::cout << " nfev: " << nfev_ << std::endl;
    std::cout << " nhev: " << nhev_ << std::endl;
  }

  /**
   * Run the optimzation algorithm for niter iterations or until the
   * stop criterion is satisfied
   */
  virtual void run() { run(maxiter_ - iter_number_); }
  /**
   * Set the initial func and gradient.  This can be used
   * to avoid one potential call
   */
  virtual void set_func_gradient(double f, Array<double> grad) {
    if (grad.size() != g_.size()) {
      throw std::invalid_argument("the gradient has the wrong size");
    }
    if (iter_number_ > 0) {
      std::cout << "warning: setting f and grad after the first iteration.This "
                   "is dangerous.\n";
    }
    // copy the function and gradient
    f_ = f;
    g_.assign(grad);
    gradient_norm_ = norm(g_) / sqrt(g_.size());
    func_initialized_ = true;
  }

  inline virtual void reset(pele::Array<double> &x0) {
    throw std::runtime_error("GradientOptimizer::reset must be overloaded");
  }

  // functions for setting the parameters
  inline void set_tol(double tol) { tol_ = tol; }
  inline void set_maxstep(double maxstep) { maxstep_ = maxstep; }
  inline void set_max_iter(int max_iter) { maxiter_ = max_iter; }
  inline void set_iprint(int iprint) { iprint_ = iprint; }
  inline void set_verbosity(int verbosity) { verbosity_ = verbosity; }
  // functions for accessing the internals for linesearches
  inline void set_f(double f_in) { f_ = f_in; }
  inline void set_rms(double rms_in) { gradient_norm_ = rms_in; }
  inline void set_potential(std::shared_ptr<pele::BasePotential> potential) {
    potential_ = potential;
  }
  inline void set_tolerance(double tol) { tol_ = tol; }

  // make sure that the allocated x is not changed externally
  // this can cause memory issues
  inline void set_x(pele::Array<double> x) { x_.assign(x); }
  // functions for accessing the status of the optimizer
  virtual inline Array<double> get_x() const { return x_; }
  inline Array<double> get_g() const { return g_; }
  inline double get_f() const { return f_; }
  inline double get_rms() const { return gradient_norm_; }
  inline int get_nfev() const { return nfev_; }
  inline int get_nhev() const { return nhev_; }
  inline int get_niter() const { return iter_number_; }
  inline int get_maxiter() const { return maxiter_; }
  inline double get_maxstep() { return maxstep_; }
  inline double get_gradient_norm_() { return gradient_norm_; }
  inline double get_step_norm_() { return step_norm_; }
  inline double get_tol() const { return tol_; }
  inline bool success() { return succeeded_; }
  inline int get_verbosity_() { return verbosity_; }
  inline bool get_last_step_failed() { return last_step_failed_; }
  inline std::shared_ptr<pele::BasePotential> get_potential() {
    return potential_;
  }

  /**
   * Return true if the termination condition is satisfied, false otherwise
   */
  virtual bool stop_criterion_satisfied() {
    if (!func_initialized_) {
      initialize_func_gradient();
    }
    succeeded_ = stop_criterion->stop_criterion_satisfied();
    return succeeded_;
  }

 protected:
  // Inactive array-gradient implementation.
  // void compute_func_gradient(Array<double> x, double & func,
  //         Array<double> gradient)
  // {
  //     nfev_ += 1;

  //     // pass the arrays to the potential
  //     func = potential_->get_energy_gradient(x, gradient);
  // }

  /** Compute the objective value and fill the supplied gradient accumulators. */
  virtual void compute_func_gradient(
      Array<double> x, double &func,
      std::vector<xsum_small_accumulator> &gradient) {
    nfev_ += 1;

    // pass the arrays to the potential
    func = potential_->get_energy_gradient(x, gradient);
  }
  /**
   * compute the initial func and gradient
   */
  virtual void initialize_func_gradient() {
    // compute the func and gradient at the current locations
    // and store them
    compute_func_gradient(x_, f_, g_);
    gradient_norm_ = norm(g_) / sqrt(x_.size());
    func_initialized_ = true;
  }

  virtual void update_costly_trajectory_info() {
    throw std::runtime_error(
        "GradientOptimizer::update_costly_trajectory_info must be overloaded");
  }
  virtual void update_cheap_trajectory_info() {
    throw std::runtime_error(
        "GradientOptimizer::update_costly_trajectory_info must be overloaded");
  }

 public:
  /**
   * functions for LineSearch
   */
  virtual void compute_func_gradient(Array<double> x, double &func,
                                     Array<double> gradient) {
    nfev_ += 1;

    // pass the arrays to the potential
    func = potential_->get_energy_gradient(x, gradient);
  }

  /**
   * Compute the norm defined by the potential
   */
  double compute_pot_norm(Array<double> x) {
    return potential_->compute_norm(x);
  }
};

/**
 * @brief Abstract base class for optimizers with a concept of time.
 */
class ODEBasedOptimizer : public GradientOptimizer {
 private:
  // Need to include others
  std::vector<double> cheap_times;  // times at which we save the cheap info
  std::vector<double> gradient_norms;
  std::vector<double> distances;
  std::vector<double> energies;

  std::vector<double> costly_times;  // times at which we save the costly info
  std::vector<std::vector<double>> coordinates;
  std::vector<std::vector<double>> gradients;

 protected:
  double time_;

  void update_costly_trajectory_info() {
    costly_times.push_back(time_);

    // we change these to std::vector<double> to
    // convert easily from the cython side
    gradients.push_back(g_.to_std_vector_copy());
    coordinates.push_back(x_.to_std_vector_copy());
  }

  void update_cheap_trajectory_info() {
    cheap_times.push_back(time_);
    gradient_norms.push_back(gradient_norm_);
    // needs to be updated by the derived class
    distances.push_back(step_norm_);
    energies.push_back(f_);
  }

 public:
  ODEBasedOptimizer(
      std::shared_ptr<pele::BasePotential> potential,
      const pele::Array<double> x0, double tol = 1e-4,
      bool save_trajectory = false, int iterations_before_save = 1,
      StopCriterionType stop_criterion = StopCriterionType::GRADIENT)
      : GradientOptimizer(potential, x0, tol, save_trajectory,
                          iterations_before_save, stop_criterion),
        time_(0) {}

  virtual ~ODEBasedOptimizer() {}

  // copy constructor raises exception that it's not implemented
  ODEBasedOptimizer(const ODEBasedOptimizer &) : GradientOptimizer() {
    throw std::runtime_error(
        "ODEBasedOptimizer copy constructor not implemented");
  }

  ODEBasedOptimizer() : GradientOptimizer() {
    throw std::runtime_error(
        "ODEBasedOptimizer default constructor not "
        "implemented. Check derived class function calls");
  }
  std::vector<double> get_time_trajectory() { return cheap_times; }
  std::vector<double> get_gradient_norm_trajectory() { return gradient_norms; }
  std::vector<double> get_distance_trajectory() { return distances; }
  std::vector<double> get_energy_trajectory() { return energies; }

  std::vector<double> get_costly_time_trajectory() { return costly_times; }
  std::vector<std::vector<double>> get_coordinate_trajectory() {
    return coordinates;
  }
  std::vector<std::vector<double>> get_gradient_trajectory() {
    return gradients;
  }
};

inline bool GradientStopCriterion::stop_criterion_satisfied() {
  if (opt_->get_gradient_norm_() < tol_) {
    return true;
  }
  return false;
}

inline bool StepnormStopCriterion::stop_criterion_satisfied() {
  if (opt_->get_gradient_norm_() < gradient_tol_) {
    if (opt_->get_step_norm_() < step_norm_tol_) {
      return true;
    }
    return false;
  }
  return false;
}

inline bool NewtonStopCriterion::stop_criterion_satisfied() {
  if (opt_->get_gradient_norm_() > gradient_tol_) {
    return false;
  }
  Array<double> hessian_pele(hessian_.data(), hessian_.size());
  Array<double> gradient_pele(gradient.data(), gradient.size());
  potential_->get_energy_gradient_hessian(opt_->get_x(), gradient_pele,
                                          hessian_pele);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es =
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(hessian_);
  eigenvalues_ = es.eigenvalues();
  eigenvectors_ = es.eigenvectors();
  double abs_min_eigenvalue = eigenvalues_.minCoeff();
  if (abs_min_eigenvalue < 0) {
    abs_min_eigenvalue = -abs_min_eigenvalue;
  } else {
    abs_min_eigenvalue = 0;
  }
  double average_eigenvalue = eigenvalues_.mean();
  double offset = std::max(offset_factor_ * std::abs(average_eigenvalue),
                           2 * abs_min_eigenvalue);
  hessian_.diagonal().array() += offset;
  newton_step_ = -hessian_.ldlt().solve(gradient);
  double newton_step_norm = newton_step_.norm();
  if (newton_step_norm > newton_tol_) {
    return false;
  }
  return true;
}
inline NewtonStopCriterion::NewtonStopCriterion(double tol,
                                                GradientOptimizer *opt)
    : gradient_tol_(tol * 1e-3),
      newton_tol_(tol),
      offset_factor_(1e-1),
      opt_(opt),
      hessian_(opt_->get_x().size(), opt_->get_x().size()),
      gradient(opt_->get_x().size()),
      eigenvalues_(opt_->get_x().size()),
      newton_step_(opt_->get_x().size()),
      eigenvectors_(opt_->get_x().size(), opt_->get_x().size()),
      potential_(opt_->get_potential()) {}

}  // namespace pele

#endif
