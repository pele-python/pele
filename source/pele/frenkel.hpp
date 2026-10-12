/******************************+
 * Modified Example from pele to create a frenkel potential
 *
 * The python bindings for this potential are in frenkel.pyx
 *
 */

#ifndef PYGMIN_FRENKEL_H
#define PYGMIN_FRENKEL_H

#include <memory>

#include "atomlist_potential.hpp"
#include "base_interaction.hpp"
#include "cell_list_potential.hpp"
#include "distance.hpp"
#include "simple_pairwise_potential.hpp"

namespace pele {

/**
 * Pairwise interaction for Frenkel potential
 * as defined in equation 3 in https://arxiv.org/abs/1910.05746
 */
struct frenkel_interaction : BaseInteraction {
  double const _eps;
  double const _sig2;
  double const _rcut2;
  double const _prefactor;
  frenkel_interaction(double sig, double eps, double rcut)
      : _eps(eps),
        _sig2(sig * sig),
        _rcut2(rcut * rcut),
        _prefactor(2 * _eps * (_rcut2 / _sig2) *
                   pow(3 / (2 * ((_rcut2 / _sig2) - 1)), 3)) {}

  /* calculate energy from distance squared */
  double inline energy(double r2, double dij) const {
    if (r2 >= _rcut2) {
      return 0.;
    }
    double ir2 = 1.0 / r2;
    double cutoff_factor = _rcut2 * ir2 - 1;
    double sigma_factor = _sig2 * ir2 - 1;
    return _prefactor * cutoff_factor * cutoff_factor * sigma_factor;
  }

  /* calculate energy and gradient from distance squared, gradient is in
   * |g|/|rij| */
  double inline energy_gradient(double r2, double *gij, double dij) const {
    double ir2 = 1.0 / r2;
    double cutoff_factor = _rcut2 * ir2 - 1;
    double sigma_factor = _sig2 * ir2 - 1;
    *gij = -_prefactor *
           (2 * _sig2 * cutoff_factor * cutoff_factor +
            4 * _rcut2 * cutoff_factor * sigma_factor) *
           ir2 * ir2;
    double energy = _prefactor * cutoff_factor * cutoff_factor * sigma_factor;
    return energy;
  }

  double inline energy_gradient_hessian(double r2, double *gij, double *hij,
                                        double dij) const {
    double ir2 = 1.0 / r2;
    double cutoff_factor = _rcut2 * ir2 - 1;
    double sigma_factor = _sig2 * ir2 - 1;
    *gij = -_prefactor *
           (2 * _sig2 * cutoff_factor * cutoff_factor +
            4 * _rcut2 * cutoff_factor * sigma_factor) *
           ir2 * ir2;
    *hij = _prefactor * (16 * _rcut2 * cutoff_factor * _sig2 * ir2 +
                         6 * cutoff_factor * cutoff_factor * _sig2 +
                         8 * _rcut2 * _rcut2 * sigma_factor * ir2 +
                         12 * _rcut2 * cutoff_factor * sigma_factor);
    double energy = _prefactor * cutoff_factor * cutoff_factor * sigma_factor;
    return energy;
  }
};

class Frenkel : public pele::SimplePairwisePotential<frenkel_interaction> {
 public:
  Frenkel(double sig, double eps, double rcut)
      : SimplePairwisePotential<frenkel_interaction>(
            std::make_shared<frenkel_interaction>(sig, eps, rcut)) {}
};

class FrenkelPeriodic : public SimplePairwisePotential<frenkel_interaction,
                                                       periodic_distance<3>> {
 public:
  FrenkelPeriodic(double sig, double eps, double rcut,
                  Array<double> const boxvec)
      : SimplePairwisePotential<frenkel_interaction, periodic_distance<3>>(
            std::make_shared<frenkel_interaction>(sig, eps, rcut),
            std::make_shared<periodic_distance<3>>(boxvec)) {}
};

/**
 * Pairwise Frenkel potential with periodic distances, using cell lists.
 */

template <size_t ndim>
class FrenkelPeriodicCellLists
    : public CellListPotential<frenkel_interaction, periodic_distance<ndim>> {
 public:
  FrenkelPeriodicCellLists(double sig, double eps, double rcut,
                           Array<double> const boxvec, double ncellx_scale)
      : CellListPotential<frenkel_interaction, periodic_distance<ndim>>(
            std::make_shared<frenkel_interaction>(sig, eps, rcut),
            std::make_shared<periodic_distance<ndim>>(boxvec), boxvec, rcut,
            ncellx_scale) {}
};

}  // namespace pele

#endif

//
//

/**
 * We now combine the components (interaction, looping method, distance
 * function) into defined classes which is the c++ potential class
 */
// class Frenkel : public pele::SimplePairwisePotential<frenkel_interaction,
// cartesian_distance<3> > { public:
//   Frenkel(double sig, double eps, double rcut)
//     : SimplePairwisePotential<frenkel_interaction, cartesian_distance<3> >
//     (std::make_shared<frenkel_interaction>(sig, eps, rcut)) {}
// };
