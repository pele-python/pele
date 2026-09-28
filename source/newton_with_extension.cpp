#include "pele/newton_with_extension.hpp"
#include <memory>

#include "pele/array.hpp"
#include "pele/eigen_interface.hpp"
#include "pele/optimizer.hpp"
// Lapack for cholesky
extern "C" {
}

using pele::Array;
using namespace std;

// value at which we assume the value is zero
#define NUMERICAL_ZERO 1e-15
namespace pele {

NewtonWithExtendedPotential::NewtonWithExtendedPotential(
    std::shared_ptr<BasePotential> potential, const pele::Array<double> &x0,
    double tol, std::shared_ptr<BasePotential> potential_extension,
    double translation_offset, double max_step)
    : GradientOptimizer(potential, x0, tol),
      _hessian(x0.size(), x0.size()),
      _gradient(x0.size()),
      _x(x0.size()),
      _x_old(x0.size()),
      _gradient_old(x0.size()),
      _line_search(this, 1.0),
      _translation_offset(translation_offset),
      _max_step(max_step),
      _potential_extension(potential_extension) {
  // Setup the extended potential
  _potential_extension = potential_extension;
  _extended_potential_wrapper =
      std::make_shared<ExtendedPotential>(potential, potential_extension);
  set_potential(
      _potential_extension);  // write pele array data into the Eigen array

  // Save coordinates as an Eigen vector
  for (size_t i = 0; i < x0.size(); ++i) {
    _x(i) = x0[i];
  }
}

void NewtonWithExtendedPotential::reset(pele::Array<double> x) {
  // write pele array data into the Eigen array
  for (size_t i = 0; i < x.size(); ++i) {
    _x(i) = x[i];
    // TODO: Hacky, find a way to remove pele arrays in the base class
    x_[i] = x[i];
  }
  initialize_func_gradient();
}

void NewtonWithExtendedPotential::one_iteration() {
  // copy in from pele to make sure initialization is taken care of
  eig_eq_pele(_x_old, x_);
  eig_eq_pele(_gradient_old, g_);

  // Wrap pele array data into the Eigen array
  Array<double> hessian_pele(_hessian.data(),
                             _hessian.rows() * _hessian.cols());
  Array<double> gradient_pele(_gradient.data(), _gradient.size());
  Array<double> x_pele(_x.data(), _x.size());

  if (hessian_calculated) {
    // do nothing
  } else {
    _energy = potential_->get_energy_gradient_hessian(x_pele, gradient_pele,
                                                      hessian_pele);
    nhev_ += 1;
  }

  // assign negative hessian;
  _hessian = -_hessian;

  double starting_norm = _step.norm();

  _x = _x_old + _step;

  // use line search to find the step size
  Array<double> gradient_old_pele(_gradient_old.data(), _gradient_old.size());
  Array<double> x_old_pele(_x_old.data(), _x_old.size());

  // Execution of the line search TODO: Clean up
  _line_search.set_xold_gold_(x_old_pele, gradient_old_pele);
  _line_search.set_g_f_ptr(gradient_pele);
  Array<double> step_pele(_step.data(), _step.size());
  double stepnorm = _line_search.line_search(x_pele, step_pele);

  if (stepnorm / starting_norm < 1e-10) {
    throw std::runtime_error(
        "rescaled step decreased by too much. Newton "
        "might be in the wrong direction");
  }
  // Careful about assignment
  x_.assign(x_pele);
  g_.assign(gradient_pele);

  iter_number_ += 1;
}

}  // namespace pele
