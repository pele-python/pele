/**
 * @file hessian_translations.hpp
 * @author Praharsh Suryadevara
 * @brief File contains functions that generate translation symmetries and
 * matrices.
 * @version 0.1
 * @date 2022-05-20
 *
 *
 */

#ifndef TRANSLATION_SYMMETRY_MATRICES_HPP
#define TRANSLATION_SYMMETRY_MATRICES_HPP

#include <Eigen/Dense>
#include <cstddef>
#include <pele/array.hpp>
#include <stdexcept>
#include <vector>

namespace pele {

inline void generate_2d_symmetries(size_t n_particles,
                                   Eigen::MatrixXd &projector_x,
                                   Eigen::MatrixXd &projector_y) {
  projector_x.setZero();
  projector_y.setZero();

  if (projector_x.rows() != 2 * n_particles) {
    throw std::logic_error("projector_x.rows() != 2 * n_particles");
  }
  if (projector_y.rows() != 2 * n_particles) {
    throw std::logic_error("projector_y.rows() != 2 * n_particles");
  }

  for (size_t i = 0; i < n_particles * 2; i++) {
    for (size_t j = i + 1; j < n_particles * 2; j++) {
      if (i % 2 == 0 && j % 2 == 0) {
        projector_x(i, j) = 1;
        projector_x(j, i) = 1;
      }
      if (i % 2 == 1 && j % 2 == 1) {
        projector_y(i, j) = 1;
        projector_y(j, i) = 1;
      }
    }
  }
}

inline void generate_3d_symmetries(size_t n_particles,
                                   Eigen::MatrixXd &projector_x,
                                   Eigen::MatrixXd &projector_y,
                                   Eigen::MatrixXd &projector_z) {
  projector_x.setZero();
  projector_y.setZero();
  projector_z.setZero();

  if (projector_x.rows() != 3 * n_particles) {
    throw std::logic_error("projector_x.rows() != 3 * n_particles");
  }
  if (projector_y.rows() != 3 * n_particles) {
    throw std::logic_error("projector_y.rows() != 3 * n_particles");
  }
  if (projector_z.rows() != 3 * n_particles) {
    throw std::logic_error("projector_z.rows() != 3 * n_particles");
  }

  for (size_t i = 0; i < n_particles * 3; i++) {
    for (size_t j = i + 1; j < n_particles * 3; j++) {
      if (i % 3 == 0 && j % 3 == 0) {
        projector_x(i, j) = 1;
        projector_x(j, i) = 1;
      }
      if (i % 3 == 1 && j % 3 == 1) {
        projector_y(i, j) = 1;
        projector_y(j, i) = 1;
      }
      if (i % 3 == 2 && j % 3 == 2) {
        projector_z(i, j) = 1;
        projector_z(j, i) = 1;
      }
    }
  }
}

inline std::vector<Eigen::MatrixXd> get_symmetries(size_t n_particles,
                                                   size_t n_d) {
  std::vector<Eigen::MatrixXd> symmetries;
  if (n_d == 2) {
    Eigen::MatrixXd projector_x, projector_y;
    projector_x = Eigen::MatrixXd::Zero(2 * n_particles, 2 * n_particles);
    projector_y = Eigen::MatrixXd::Zero(2 * n_particles, 2 * n_particles);
    generate_2d_symmetries(n_particles, projector_x, projector_y);
    symmetries.push_back(projector_x);
    symmetries.push_back(projector_y);
  } else if (n_d == 3) {
    Eigen::MatrixXd projector_x, projector_y, projector_z;
    projector_x = Eigen::MatrixXd::Zero(2 * n_particles, 2 * n_particles);
    projector_y = Eigen::MatrixXd::Zero(2 * n_particles, 2 * n_particles);
    generate_3d_symmetries(n_particles, projector_x, projector_y, projector_z);
    symmetries.push_back(projector_x);
    symmetries.push_back(projector_y);
    symmetries.push_back(projector_z);
  } else {
    throw std::logic_error("n_d must be 2 or 3");
  }
  return symmetries;
}

/**
 * @brief  adds a translation offset to the hessian. This should take care of
 * translational symmetries
 * @param  hessian the hessian to be modified
 * @param  offset the offset to be added
 */
inline void add_translation_offset_2d(Eigen::MatrixXd &hessian, double offset) {
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
  hessian.diagonal().array() += offset;
}

/**
 * @brief  adds a translation offset to the hessian. This should take care of
 * translational symmetries. offset[i * hessian.rows() + j] is the offset at the
 * i, jth element, only the upper triangle is added since
 * the hessian is symmetric
 * @param  hessian the hessian to be modified
 * @param  offset the offset to be added this is a pele array in this case
 */
inline void add_symmetry_offset(Eigen::MatrixXd &hessian,
                                const Array<double> &offset) {
  for (size_t i = 0; i < hessian.rows(); ++i) {
    for (size_t j = 0; j < hessian.cols(); ++j) {
      hessian(i, j) += offset[i * hessian.rows() + j];
    }
  }
}

}  // namespace pele
#endif  // TRANSLATION_SYMMETRY_MATRICES_HPP
