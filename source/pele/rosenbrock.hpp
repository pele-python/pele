#ifndef _PELE_TEST_FUNCS_H
#define _PELE_TEST_FUNCS_H

#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <numbers>
#include <numeric>
#include <vector>

#include "array.hpp"
#include "base_potential.hpp"

namespace pele {
/**
 * RosenBrock function with defaults for testing purposes
 */
class RosenBrock : public BasePotential {
 private:
  double m_a;
  double m_b;

 public:
  RosenBrock(double a = 1, double b = 100) : m_a(a), m_b(b){};
  virtual ~RosenBrock(){};
  inline double get_energy(Array<double> const &x) {
    return (m_a - x[0]) * (m_a - x[0]) +
           m_b * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]);
  }
  inline double get_energy_gradient(Array<double> const &x,
                                    Array<double> &grad) {
    grad.assign(0);
    grad[0] = 4 * m_b * x[0] * x[0] * x[0] - 4 * m_b * x[0] * x[1] +
              2 * m_a * x[0] - 2 * m_a;
    grad[1] = 2 * m_b * (x[1] - x[0] * x[0]);
    return (m_a - x[0]) * (m_a - x[0]) +
           m_b * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]);
  };
};

/**
 * Saddle point function for optimizer testing purposes
 */
class Saddle : public BasePotential {
 private:
  double m_a2;
  double m_b2;
  double m_a4;
  double m_b4;

 public:
  Saddle(double a2 = 2, double b2 = 2, double a4 = 1, double b4 = 1)
      : m_a2(a2), m_b2(b2), m_a4(a4), m_b4(b4){};
  virtual ~Saddle(){};
  inline double get_energy(Array<double> const &x) {
    return m_a2 * x[0] * x[0] - m_b2 * x[1] * x[1] +
           m_a4 * x[0] * x[0] * x[0] * x[0] + m_a4 * x[1] * x[1] * x[1] * x[1];
  }
  inline double get_energy_gradient(Array<double> const &x,
                                    Array<double> &grad) {
    grad.assign(0);
    grad[0] = 2 * m_a2 * x[0] + 4 * m_b4 * x[0] * x[0] * x[0];
    grad[1] = -2 * m_b2 * x[1] + 4 * m_b4 * x[1] * x[1] * x[1];
    return m_a2 * x[0] * x[0] - m_b2 * x[1] * x[1] +
           m_a4 * x[0] * x[0] * x[0] * x[0] + m_a4 * x[1] * x[1] * x[1] * x[1];
  };
};

/**
 * 1d x cubed function for testing optimizers. WARNING stop minimization before
 * the energy goes too low.
 */
class XCube : public BasePotential {
 public:
  XCube(){};
  virtual ~XCube(){};
  inline double get_energy(Array<double> const &x) {
    return x[0] * x[0] * x[0];
  }
  inline double get_energy_gradient(Array<double> const &x,
                                    Array<double> &grad) {
    grad.assign(0);
    grad[0] = 3 * x[0] * x[0];
    std::cout << grad[0] << " graaadient \n";
    return x[0] * x[0] * x[0];
  };
};

/**
 * Flat Harmonic Potential. Flat in one dimension so that the hessian is
 * singular. Useful for testing optimizers where the singularity is handled.
 * The Potential is flat in the x+y direction. and the origin is at 0, 0.
 */

class FlatHarmonic : public BasePotential {
 private:
  double m_a;

 public:
  FlatHarmonic(double a = 1) : m_a(a){};
  virtual ~FlatHarmonic(){};
  inline double get_energy(Array<double> const &x) { return m_a * x[0] * x[0]; }
  inline double get_energy_gradient(Array<double> const &x,
                                    Array<double> &grad) {
    grad.assign(0);
    grad[0] = 2 * m_a * (x[0] - x[1]);
    grad[1] = -grad[0];
    return m_a * x[0] * x[0];
  };
  inline double get_energy_gradient_hessian(Array<double> const &x,
                                            Array<double> &grad,
                                            Array<double> &hess) {
    grad.assign(0);
    hess.assign(0);
    grad[0] = 2 * m_a * (x[0] - x[1]);
    grad[1] = -grad[0];
    hess[0] = 2 * m_a;
    hess[1] = -hess[0];
    hess[2] = -hess[0];
    hess[3] = hess[0];
    return m_a * x[0] * x[0];
  };
};
/**
 * Negative cosine sum function for testing basin volume calculations.
 * Yields a hypercubic basin with volume period^dim.
 * We also add a power to the energy to ensure that the hessian isn't diagonal.
 * The full potential is given by
 * \f[
 * V(x) = \left(\mathrm{offset} + \sum_i
 * \left[1-\cos(2\pi x_i/\mathrm{period})\right]\right)^{\mathrm{power}}
 * \f]
 * The offset ensures that the hessian is smooth if power goes below 0.
 */
class PoweredCosineSum : public BasePotential {
 private:
  double _period_factor;
  // precompute cos and sin values
  Array<double> _cos_values;
  Array<double> _sin_values;
  double _power;
  double _offset;
  void _precompute_cos_sin(Array<double> const &x) {
    for (size_t i = 0; i < x.size(); ++i) {
      _cos_values[i] = std::cos(_period_factor * x[i]);
      _sin_values[i] = std::sin(_period_factor * x[i]);
    }
  }

 public:
  PoweredCosineSum(size_t dim, double period = 1, double pow = 0.5,
                   double offset = 1.0)
      : _period_factor(2 * std::numbers::pi / period),
        _cos_values(dim, 0.0),
        _sin_values(dim, 0.0),
        _power(pow),
        _offset(offset){};
  virtual ~PoweredCosineSum(){};
  inline double get_energy(Array<double> const &x) {
    if (_cos_values.size() != x.size()) {
      throw std::runtime_error(
          "dim passed to cosine sume not the same size as input");
    }
    _precompute_cos_sin(x);
    double f_x = x.size() + _offset;
    f_x -= std::accumulate(_cos_values.begin(), _cos_values.end(), 0.0);
    return std::pow(f_x, _power);
  }

  double add_energy_gradient(const Array<double> &x, Array<double> &grad) {
    if (_cos_values.size() != x.size()) {
      throw std::runtime_error(
          "dim passed to cosine sume not the same size as input");
    }
    _precompute_cos_sin(x);
    double f_x = x.size() + _offset;
    f_x -= std::accumulate(_cos_values.begin(), _cos_values.end(), 0.0);
    // m for minus
    double power_m_1_term = std::pow(f_x, _power - 1);

    for (size_t i = 0; i < x.size(); ++i) {
      grad[i] += _power * power_m_1_term * _period_factor * _sin_values[i];
    }
    return power_m_1_term * f_x;
  }

  double get_energy_gradient(const Array<double> &x, Array<double> &grad) {
    grad.assign(0);
    return add_energy_gradient(x, grad);
  }

  void add_hessian(const Array<double> &x, Array<double> &hess) {
    if (_cos_values.size() != x.size()) {
      throw std::runtime_error(
          "dim passed to cosine sume not the same size as input");
    }
    _precompute_cos_sin(x);
    double f_x = x.size() + _offset;
    f_x -= std::accumulate(_cos_values.begin(), _cos_values.end(), 0.0);
    // m for minus
    double power_m_2 = std::pow(f_x, _power - 2);
    for (size_t i = 0; i < x.size(); ++i) {
      for (size_t j = 0; j < x.size(); ++j) {
        double common_term = _power * (_power - 1) * power_m_2 *
                             _period_factor * _period_factor * _sin_values[i] *
                             _sin_values[j];
        if (i == j) {
          hess[i * x.size() + j] +=
              common_term + _power * power_m_2 * f_x * _period_factor *
                                _period_factor * _cos_values[i];
        } else {
          hess[i * x.size() + j] += common_term;
        }
      }
    }
  }

  inline double add_energy_gradient_hessian(Array<double> const &x,
                                            Array<double> &grad,
                                            Array<double> &hess) {
    if (_cos_values.size() != x.size()) {
      throw std::runtime_error(
          "dim passed to cosine sume not the same size as input");
    }
    _precompute_cos_sin(x);
    double f_x = x.size() + _offset;
    f_x -= std::accumulate(_cos_values.begin(), _cos_values.end(), 0.0);
    // m for minus
    double power_m_2 = std::pow(f_x, _power - 2);
    for (size_t i = 0; i < x.size(); ++i) {
      grad[i] += _power * power_m_2 * f_x * _period_factor * _sin_values[i];
    }

    for (size_t i = 0; i < x.size(); ++i) {
      for (size_t j = 0; j < x.size(); ++j) {
        double common_term = _power * (_power - 1) * power_m_2 *
                             _period_factor * _period_factor * _sin_values[i] *
                             _sin_values[j];
        if (i == j) {
          hess[i * x.size() + j] +=
              common_term + _power * power_m_2 * f_x * _period_factor *
                                _period_factor * _cos_values[i];
        } else {
          hess[i * x.size() + j] += common_term;
        }
      }
    }
    return power_m_2 * f_x * f_x;
  }

  double get_energy_gradient_hessian(Array<double> const &x,
                                     Array<double> &grad, Array<double> &hess) {
    grad.assign(0);
    hess.assign(0);
    return add_energy_gradient_hessian(x, grad, hess);
  }
};

}  // namespace pele
#endif
