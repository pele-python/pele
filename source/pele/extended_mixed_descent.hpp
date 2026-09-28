#ifndef _PELE_EXTENDED_MIXED_DESCENT_H__
#define _PELE_EXTENDED_MIXED_DESCENT_H__

#include "array.hpp"
#include "base_potential.hpp"
#include "preprocessor_directives.hpp"

// #define EIGEN_USE_MKL_ALL
// Eigen linear algebra library
#include <Eigen/Dense>
#include <array>
#include <fstream>
#include <memory>
#include <sundials/sundials_context.h>

#include "eigen_interface.hpp"
#include "pele/combine_potentials.hpp"
#include "pele/lbfgs.hpp"

// Lapack for cholesky
extern "C" {
#include <lapack.h>
}

// line search methods
#include <cvode/cvode.h>               /* access to CVODE                 */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector       */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix       */

#include "backtracking.hpp"
#include "bracketing.hpp"
#include "cvode.hpp"
#include "linesearch.hpp"
#include "more_thuente.hpp"
#include "nwpele.hpp"
#include "optimizer.hpp"

extern "C" {
#include "xsum.h"
}

namespace pele {

class ExtendedMixedOptimizer : public GradientOptimizer {
 private:
  /**
   * H0 is the initial estimate for the inverse hessian. This is as small as
   * possible, for making a good inverse hessian estimate next time.
   */
  /**
   * User_Data object. This is important for mxopt because it stores the number
   * of getEnergyGradient and getEnergyGradientHessian evaluations. We're saving
   * them here instead of the optimizer because we'd need to keep a double count
   * since SUNDIALS needs this object.
   */
  UserData_ udata;
  void *cvode_mem; /* CVODE memory         */
  size_t N_size;
  SUNMatrix A;
  SUNLinearSolver LS;
  SUNContext sunctx;
  double t0;
  double tN;
  N_Vector x0_N;
  double rtol;
  double atol;

  // return val from CVODE
  int retval;
  // instatiate if you want to print to file
#if PRINT_TO_FILE == 1
  std::ofstream trajectory_file;
  std::ofstream hessian_eigvals_file;
#endif
  /**
   * extension of the original potential that ensures the hessian is positive
   * definite. This has a switch to switch on/off the extension.
   */
  std::shared_ptr<pele::ExtendedPotential> extended_potential;

  Array<double> xold;  //!< Coordinates before taking a step
  Array<double> gold;  //!< Gradient before taking a step
  Array<double> step;  //!< Step size and direction

  /**
   * Coordinate for phase 1. Helps us revert back if newton fails.
   */
  Array<double> xold_old;  //!< Save for backtracking if newton fails. TODO:
                           //!< think of better name

  /**
   * @brief Last CVODE position to return to if newton fails.
   */
  Array<double> x_last_cvode;

  double inv_sqrt_size;  //!< The inverse square root the the number of
                         //!< components
  // Preconditioning
  int T_;  // number of steps after which the lowest eigenvalues are
           // recalculated in the first phase
  // std::shared_ptr<Eigen::ColPivHouseholderQR<Eigen::MatrixXd>> solver;
  // solver for H x = b
  Eigen::VectorXd update_solver(
      Eigen::VectorXd r);  // updates the solver with the new hessian
  char uplo;               /* We ask LAPACK for the lower diagonal matrix L */
  int info;                /* "Info" return value, used for error-checking */

  // Calculates hess + delta I where delta makes the new eigenvalue positive
  Eigen::MatrixXd hessian;
  Eigen::MatrixXd hessian_copy_for_cholesky;

  double *hess_data;
  bool use_phase_1;
  bool phase_2_failed_;  // if true, we've failed in phase 2 and need to revert
                         // back to where we were in phase 1
  // what phase was used in previous step
  bool prev_phase_is_phase1;
  /**
   * number of phase 1 steps
   */
  size_t n_phase_1_steps;
  /**
   * number of phase 2 steps
   */
  size_t n_phase_2_steps;

  /**
   * @brief Failed phase 2 steps
   */
  size_t n_failed_phase_2_steps;

  /**
   * tolerance for convexity. the smaller, the more convex the problem
   * needs to be before switching to newton
   */
  double conv_tol_;
  // Need to refactor line searches
  BacktrackingLineSearch line_search_method;

  bool iterative_;  // if true, use iterative solver for hessian solve
  HessianType hessian_type_;

  const Array<double>
      m_global_symmetry_offset;  // global symmetry offset for hessian
                                 // to account for numerical issues with
                                 // cholesky finding negative eigenvalues
                                 // If size zero defaults to using 2d
                                 // translational symmetries to maintain
                                 // backwards compatibility TODO: deprecate
                                 // this default

 public:
  /**
   * Constructor
   */
  ExtendedMixedOptimizer(
      std::shared_ptr<pele::BasePotential> potential,
      const pele::Array<double> x0,
      std::shared_ptr<pele::BasePotential> potential_extension = nullptr,
      double tol = 1e-4, int T = 10, double step = 1, double conv_tol = 1e-8,
      double rtol = 1e-5, double atol = 1e-5, bool iterative = false,
      Array<double> global_symmetry_offset = Array<double>(0),
      bool save_trajectory = false, int iterations_before_save = 1,
      StopCriterionType stop_criterion = StopCriterionType::GRADIENT);
  /**
   * Destructor
   */
  virtual ~ExtendedMixedOptimizer();

  void free_cvode_objects();

  /**
   * Do one iteration iteration of the optimization algorithm
   */
  void one_iteration();

  // functions for accessing the results

  /**
   * reset the optimizer to start a new minimization from x0
   */
  void reset(pele::Array<double> &x0);

  /**
   * @brief reset cvode data
   *
   */
  void reset_cvode();

  /**
   * @brief reset the newton method
   *
   */
  void reset_newton();

  inline int get_nhev() const { return udata.nhev; }
  inline int get_nhev_extended() const {
    return extended_potential->get_nhev_extension();
  }
  inline int get_n_phase_1_steps() { return n_phase_1_steps; }
  inline int get_n_phase_2_steps() { return n_phase_2_steps; }
  inline int get_n_failed_phase_2_steps() { return n_failed_phase_2_steps; }

 private:
  bool hessian_calculated;  // checks whether the hessian has been calculated
                            // for updating.
  void compute_phase_1_step();
  void compute_phase_2_step();
  bool convexity_check();
  void get_hess(Eigen::MatrixXd &hess);
  void get_hess_extended(Eigen::MatrixXd &hess);
  void setup_cvode();

  void update_H0_(Array<double> x_old, Array<double> &g_old,
                  Array<double> x_new, Array<double> &g_new);
};

}  // namespace pele

#endif  // _PELE_EXTENDED_MIXED_DESCENT_H__
