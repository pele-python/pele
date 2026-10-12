#ifndef _PELE_LINESEARCH_H__
#define _PELE_LINESEARCH_H__

/**
 * Should include line search methods for use in various optimizers
 * refer to Nocedal and Wright Chapter 3 for information.
 * preferably should be header only,
 *
 *
 *
 */

#include <memory>
#include <vector>

#include "array.hpp"
#include "base_potential.hpp"
#include "optimizer.hpp"

namespace pele {

/**
 * Abstract base class for Line Search Methods in optimizers. All Line search
 * methods should ideally derive from this class
 */
class LineSearch {
 protected:
  // // stored energy and gradient at the computed stepsize
  Array<double> end_gradient;
  double end_energy;
  double stepnorm;
  // Base potential pointer
  double max_stepnorm_;
  // Optimizer class
  GradientOptimizer *opt_;

 public:
  LineSearch(GradientOptimizer *opt, double max_stepnorm)
      : opt_(opt), max_stepnorm_(max_stepnorm){};
  virtual ~LineSearch(){};

  /** Compute the step size given a step and a position. */
  virtual double line_search(Array<double> &x, Array<double> step) = 0;

  // // variable settings
  inline void set_max_stepnorm_(double max_stepnorm_in) {
    max_stepnorm_ = max_stepnorm_in;
  }
};

/**
 * line search method as implemented in backtracking line search in pele
 * originally for LBFGS: as defined by Daniel Asenjo prevents the energy from
 * rising above a certain value as compared to the usual backtracking linesearch
 * method which implements sufficient decrease condition May also fail if the
 * convergence is quadratic
 *
 * ALso it is assumed that any potential methods used by this would be defined
 * from an interface to the optimizer so that the evaluations can be tracked
 */

class OldLineSearch : public LineSearch {
 public:
  OldLineSearch(GradientOptimizer *opt, double maxstepnorm,
                double max_f_rise = 1e-4, double nred_max = 10,
                bool use_relative_f = true)
      : LineSearch(opt, maxstepnorm),
        max_f_rise_(max_f_rise),
        nred_max_(nred_max),
        use_relative_f_(use_relative_f),
        verbosity_(opt->get_verbosity_()){

        };

 private:
  // maximum the function can rise to
  double max_f_rise_;
  // maximum number of backtracking iterations
  int nred_max_;
  // whether to use relative f or absolute f in determining how much larger the
  // function can go
  bool use_relative_f_;
  Array<double> xold;
  // verbosity, should be inherited from the optimizer
  int verbosity_;
  Array<double> x_;
  Array<double> g_;
  Array<double> gold_;
  Array<double> xold_;
  // Set an f pointer
  double *f_ptr_;
  double f_;
  // GradientOptimizer * gopt_ ; // hacky but I don't know any good way to
  // specify opt is gradient optimizer This ends up giving me a version of this
  // with version

 public:
  // I'm putting this here instead of putting it elsewhere since I ended up
  // getting an undefined symbol error which means there are some linker issues?
  double line_search(Array<double> &x, Array<double> step);
  // This is slightly convoluted because we want to set the variables without
  // setting them as args in the linesearch method.

  // helper function to assign xold and gold.

  // These should reflect the xold in the old version
  inline void set_xold_gold_(Array<double> &xold, Array<double> &gold) {
    gold_ = gold;
    xold_ = xold;
  };
  // helper function to assign g_ and f
  inline void set_g_f_ptr(Array<double> &g) { g_ = g; };
};

}  // namespace pele

#endif /* End line search methods */
