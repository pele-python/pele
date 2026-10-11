/**
 * @file generic_mixed_descent.cpp
 * @author Praharsh Suryadevara (praharsharmm@gmail.com)
 * @brief Implementations for generic_mixed_descent.hpp
 * @version 0.1
 * @date 2022-06-07
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "pele/generic_mixed_descent.hpp"
#include <lapack.h>
#include <ostream>
#include <pele/array.hpp>

using namespace std;
namespace pele {
bool GenericMixedDescent::x_in_convex_region() {
  x_.assign(opt_non_convex_->get_x());
  // Calculate hessian at x. Current version is hessian based.
  // Ideally you should be able to get the pointer to the hessian from the
  // optimizer.
  get_hess(hessian_);

  int N_int = x_.size();
  add_translation_offset_2d(hessian_, translation_offset_);

  hessian_copy_for_cholesky_ = hessian_;

  // pointer to t
  double *hessian_copy_for_cholesky_data;
  hessian_copy_for_cholesky_data = hessian_copy_for_cholesky_.data();
  int info = 0;
  dpotrf_(&uplo, &N_int, hessian_copy_for_cholesky_data, &N_int, &info
#ifdef LAPACK_FORTRAN_STRLEN_END  // 1 is the length of the uplo string,
                                  // necessary for current fortran
          ,
          1
#endif
  );

  if (info == 0) {
    // if positive definite, we are in the convex region.
    return true;
  }
  return false;
}

void GenericMixedDescent::get_hess(Eigen::MatrixXd &hessian) {
  // Does not allocate memory for hessian just wraps around the data
  Array<double> hessian_pele = Array<double>(hessian.data(), hessian.size());
  if (pot_extension_provided_) {
    extended_potential_->get_hessian(x_, hessian_pele);
  } else {
    potential_->get_hessian(x_, hessian_pele);
  }
  hessian_evaluations_ += 1;
}

void GenericMixedDescent::one_iteration() {
  if ((iter_number_ % steps_before_convex_check_ == 0 and iter_number_ > 0) and
      !in_convex_region_) {
    in_convex_region_ = x_in_convex_region();

    // if in convex region and non convex method is being used
    // switch to convex method
    if (in_convex_region_ and use_non_convex_method_) {
      use_non_convex_method_ = false;
      last_non_convex_x_.assign(x_);
      opt_convex_->set_x(x_);
    }
  }
  // take iteration
  if (use_non_convex_method_) {
    opt_non_convex_->one_iteration();
    n_phase_1_steps += 1;
    gradient_norm_ = opt_non_convex_->get_rms();

    // Debugging
#if OPTIMIZER_DEBUG_LEVEL > 2
    std::cout << "non convex step" << std::endl;
    std::cout << "rms" << gradient_norm_ << std::endl;
#endif

  } else {
    opt_convex_->one_iteration();
    gradient_norm_ = opt_convex_->get_rms();

    // Debugging
#if OPTIMIZER_DEBUG_LEVEL > 2
    std::cout << "convex step" << std::endl;
    std::cout << "rms" << gradient_norm_ << std::endl;
#endif
    if (opt_convex_->get_last_step_failed()) {
      use_non_convex_method_ = true;
      opt_non_convex_->set_x(last_non_convex_x_);
    }
  }
  // update all parameters from the optimizer
  nfev_ = opt_non_convex_->get_nfev() + opt_convex_->get_nfev();
  nhev_ = opt_non_convex_->get_nhev() + opt_convex_->get_nhev() +
          hessian_evaluations_;
  iter_number_ += 1;
}

/**
 * @brief  adds a translation offset to the hessian. This should take care of
 * translational symmetries
 * @param  hessian the hessian to be modified
 * @param  offset the offset to be added
 */
void GenericMixedDescent::add_translation_offset_2d(Eigen::MatrixXd &hessian,
                                                    double offset) {
  // factor so that the added translation operators are unitary
  double factor = 2.0 / hessian.rows();
  offset = offset * factor;

  for (size_t i = 0; i < hessian.rows(); ++i) {
    for (size_t j = i + 1; j < hessian.cols(); ++j) {
      if ((i % 2 == 0 && j % 2 == 0) || (i % 2 == 1 && j % 2 == 1)) {
        hessian(i, j) += offset;
        hessian(j, i) += offset;
      }
    }
  }
  // add the diagonal elements that were missed by the previous loop
  hessian.diagonal().array() += offset;
}
}  // namespace pele
