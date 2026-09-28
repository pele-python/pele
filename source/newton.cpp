#include "pele/newton.hpp"
#include <fstream>

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

Newton::Newton(std::shared_ptr<BasePotential> potential,
               const pele::Array<double> &x0, double tol, double threshold,
               bool use_rattler_mask)
    : GradientOptimizer(potential, x0, tol),
      _hessian(x0.size(), x0.size()),
      _gradient(x0.size()),
      _x(x0.size()),
      _x_old(x0.size()),
      _gradient_old(x0.size()),
      _line_search(this, 1.0),
      _tolerance(tol),
      _threshold(threshold),
      _nhev(0),
      _use_rattler_mask(use_rattler_mask),
      _not_rattlers(x0.size(), true),
      _rattlers_found(false)  // use step size of 1.0 for newton
{
  // write pele array data into the Eigen array
  for (size_t i = 0; i < x0.size(); ++i) {
    _x(i) = x0[i];
  }
  // // Find rattlers is not necessary here in the way we're using it.
  // if (_use_rattler_mask) {
  //   potential->find_rattlers(x0, _not_rattlers, _jammed);
  // }
}

void Newton::reset(pele::Array<double> x) {
  if (_use_rattler_mask) {
    potential_->find_rattlers(x, _not_rattlers, _jammed);
    _rattlers_found = true;
  }

  // write pele array data into the Eigen array
  for (size_t i = 0; i < x.size(); ++i) {
    _x(i) = x[i];
    // TODO: Hacky, find a way to remove pele arrays in the base class
    x_[i] = x[i];
  }
  initialize_func_gradient();
}

void Newton::one_iteration() {
  // copy in from pele to make sure initialization is taken care of
  eig_eq_pele(_x_old, x_);
  eig_eq_pele(_gradient_old, g_);

  // Wrap pele array data into the Eigen array
  Array<double> hessian_pele(_hessian.data(),
                             _hessian.rows() * _hessian.cols());
  Array<double> gradient_pele(_gradient.data(), _gradient.size());
  Array<double> x_pele(_x.data(), _x.size());

  if (_use_rattler_mask) {
    _energy = potential_->get_energy_gradient_hessian(x_pele, gradient_pele,
                                                      hessian_pele);

    _energy = potential_->get_energy_gradient_hessian_rattlers(
        x_pele, gradient_pele, hessian_pele, _not_rattlers);
  } else {
    _energy = potential_->get_energy_gradient_hessian(x_pele, gradient_pele,
                                                      hessian_pele);
  }
  _nhev += 1;

  cout << "gradient pele" << gradient_pele << endl;

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es =
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(_hessian);
  Eigen::VectorXd eigenvalues = es.eigenvalues();
  Eigen::MatrixXd eigenvectors = es.eigenvectors();
  // eigenvector corresponding to smallest eigenvalue
  Eigen::VectorXd smallest_eigenvector = eigenvectors.col(0);

  Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(_hessian);

  // replace eigenvalues with their absolute values
  for (size_t i = 0; i < eigenvalues.size(); ++i) {
    eigenvalues(i) = std::abs(eigenvalues(i));
  }

  // cout << "eigenvalues" << eigenvalues << endl;
  // postprocess inverse Eigen
  Eigen::VectorXd inv_eigenvalues(eigenvalues.size());
  for (size_t i = 0; i < eigenvalues.size(); ++i) {
    double mod_eig = abs(eigenvalues(i));
    if (mod_eig < NUMERICAL_ZERO) {
      // case for translational symmetries
      inv_eigenvalues(i) = 0.0;
    } else if (mod_eig < _threshold) {
      // case for really small eigenvalues
      inv_eigenvalues(i) = 1.0 / _threshold;
    } else {
      // case for normal eigenvalues
      inv_eigenvalues(i) = 1.0 / mod_eig;
    }
  }

  // find the largest eigenvalue
  _step = -eigenvectors *
          ((eigenvectors.transpose() * _gradient).cwiseProduct(inv_eigenvalues))
              .matrix();

  // cout << "gradient \n" << _gradient << endl;

  // cout << "step size: " << _step.norm() << endl;
  double starting_norm = _step.norm();

  // _step = smallest_eigenvector;
  // calculate the newton step
  // _step = -cod.solve(_gradient);

  _x = _x_old + _step;

  // use line search to find the step size
  Array<double> gradient_old_pele(_gradient_old.data(), _gradient_old.size());
  Array<double> x_old_pele(_x_old.data(), _x_old.size());

  _line_search.set_xold_gold_(x_old_pele, gradient_old_pele);
  _line_search.set_g_f_ptr(gradient_pele);
  Array<double> step_pele(_step.data(), _step.size());
  double stepnorm = _line_search.line_search(x_pele, step_pele);

  // cout << "step \n" << _step << endl;
  // cout << "gradient \n" << _gradient << endl;
  // cout << "hessian \n" << _hessian << endl;
  // cout << "eigenvectors \n" << eigenvectors << endl;
  // cout << "eigenvalues \n" << eigenvalues << endl;

  // cout << "starting norm" << starting_norm << endl;
  // cout << "step norm" << stepnorm << endl;
  if (stepnorm / starting_norm < 1e-10) {
    throw std::runtime_error(
        "rescaled step decreased by too much. Newton "
        "might be in the wrong direction");
  }
  cout << "starting step size" << starting_norm << endl;
  cout << "step size rescaled" << stepnorm << endl;
  // Hacky since we're not wrapping the original pele array. but we can change
  // this if this if it becomes a speed constraint/causes maintainability issues
  x_ = x_pele;
  g_ = gradient_pele;

  iter_number_ += 1;
}

void Newton::compute_func_gradient(Array<double> x, double &func,
                                   Array<double> gradient) {
  nfev_ += 1;

  // pass the arrays to the potential
  if (_use_rattler_mask) {
    func = potential_->get_energy_gradient(x, gradient, _not_rattlers);
  } else {
    func = potential_->get_energy_gradient(x, gradient);
  }
}

}  // namespace pele