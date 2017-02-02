#ifndef _PELE_CELL_LIST_POTENTIAL_H
#define _PELE_CELL_LIST_POTENTIAL_H

#include <iostream>
#include <memory>
#include <vector>
#include <utility>
#include <stdexcept>

#include "pairwise_potential_interface.h"
#include "array.h"
#include "distance.h"
#include "cell_lists.h"
#include "vecn.h"

namespace pele{

/**
 * class which accumulates the energy one pair interaction at a time
 */
template <typename pairwise_interaction, typename distance_policy>
class EnergyAccumulator {
    const static size_t m_ndim = distance_policy::_ndim;
    std::shared_ptr<pairwise_interaction> & m_interaction;
    std::shared_ptr<distance_policy> & m_dist;

public:
    double m_energy;

    EnergyAccumulator(std::shared_ptr<pairwise_interaction> & interaction,
            std::shared_ptr<distance_policy> & dist)
        : m_interaction(interaction),
          m_dist(dist),
          m_energy(0.)
    {}

    void insert_atom_pair(Array<double> const & coords, const size_t atom_i, const size_t atom_j)
    {
        double dr[m_ndim];
        m_dist->get_rij(dr, coords.data() + m_ndim * atom_i, coords.data() + m_ndim * atom_j);
        double r2 = 0;
        for (size_t k = 0; k < m_ndim; ++k) {
            r2 += dr[k] * dr[k];
        }
        m_energy += m_interaction->energy(r2, atom_i, atom_j);
    }
};

/**
 * class which accumulates the energy and gradient one pair interaction at a time
 */
template <typename pairwise_interaction, typename distance_policy>
class EnergyGradientAccumulator {
    const static size_t m_ndim = distance_policy::_ndim;
    std::shared_ptr<pairwise_interaction> & m_interaction;
    std::shared_ptr<distance_policy> & m_dist;

public:
    double m_energy;
    pele::Array<double> & m_gradient;

    EnergyGradientAccumulator(std::shared_ptr<pairwise_interaction> & interaction,
            std::shared_ptr<distance_policy> & dist, pele::Array<double> & gradient)
        : m_interaction(interaction),
          m_dist(dist),
          m_energy(0.),
          m_gradient(gradient)
    {}

    void insert_atom_pair(Array<double> const & coords, const size_t atom_i, const size_t atom_j)
    {
        const size_t xi_off = m_ndim * atom_i;
        const size_t xj_off = m_ndim * atom_j;
        double dr[m_ndim];
        m_dist->get_rij(dr, coords.data() + m_ndim * atom_i, coords.data() + m_ndim * atom_j);
        double r2 = 0;
        for (size_t k = 0; k < m_ndim; ++k) {
            r2 += dr[k] * dr[k];
        }
        double gij;
        m_energy += m_interaction->energy_gradient(r2, &gij, atom_i, atom_j);
        for (size_t k = 0; k < m_ndim; ++k) {
            m_gradient[xi_off + k] -= gij * dr[k];
        }
        for (size_t k = 0; k < m_ndim; ++k) {
            m_gradient[xj_off + k] += gij * dr[k];
        }
    }
};

/**
 * class which accumulates the energy, gradient, and Hessian one pair interaction at a time
 */
template <typename pairwise_interaction, typename distance_policy>
class EnergyGradientHessianAccumulator {
    const static size_t m_ndim = distance_policy::_ndim;
    std::shared_ptr<pairwise_interaction> & m_interaction;
    std::shared_ptr<distance_policy> & m_dist;

public:
    double m_energy;
    pele::Array<double> & m_gradient;
    pele::Array<double> & m_hessian;

    EnergyGradientHessianAccumulator(std::shared_ptr<pairwise_interaction> & interaction,
            std::shared_ptr<distance_policy> & dist, pele::Array<double> & gradient,
            pele::Array<double> & hessian)
        : m_interaction(interaction),
          m_dist(dist),
          m_energy(0.),
          m_gradient(gradient),
          m_hessian(hessian)
    {}

    void insert_atom_pair(Array<double> const & coords, const size_t atom_i, const size_t atom_j)
    {
        const size_t xi_off = m_ndim * atom_i;
        const size_t xj_off = m_ndim * atom_j;
        double dr[m_ndim];
        m_dist->get_rij(dr, coords.data() + m_ndim * atom_i, coords.data() + m_ndim * atom_j);
        double r2 = 0;
        for (size_t k = 0; k < m_ndim; ++k) {
            r2 += dr[k] * dr[k];
        }
        double gij, hij;
        m_energy += m_interaction->energy_gradient_hessian(r2, &gij, &hij, atom_i, atom_j);
        for (size_t k = 0; k < m_ndim; ++k) {
            m_gradient[xi_off + k] -= gij * dr[k];
        }
        for (size_t k = 0; k < m_ndim; ++k) {
            m_gradient[xj_off + k] += gij * dr[k];
        }
        //this part is copied from simple_pairwise_potential.h
        //(even more so than the rest)
        const size_t N = m_gradient.size();
        const size_t i1 = xi_off;
        const size_t j1 = xj_off;
        for (size_t k = 0; k < m_ndim; ++k) {
            //diagonal block - diagonal terms
            const double Hii_diag = (hij + gij) * dr[k] * dr[k] / r2 - gij;
            m_hessian[N * (i1 + k) + i1 + k] += Hii_diag;
            m_hessian[N * (j1 + k) + j1 + k] += Hii_diag;
            //off diagonal block - diagonal terms
            const double Hij_diag = -Hii_diag;
            m_hessian[N * (i1 + k) + j1 + k] = Hij_diag;
            m_hessian[N * (j1 + k) + i1 + k] = Hij_diag;
            for (size_t l = k + 1; l < m_ndim; ++l) {
                //diagonal block - off diagonal terms
                const double Hii_off = (hij + gij) * dr[k] * dr[l] / r2;
                m_hessian[N * (i1 + k) + i1 + l] += Hii_off;
                m_hessian[N * (i1 + l) + i1 + k] += Hii_off;
                m_hessian[N * (j1 + k) + j1 + l] += Hii_off;
                m_hessian[N * (j1 + l) + j1 + k] += Hii_off;
                //off diagonal block - off diagonal terms
                const double Hij_off = -Hii_off;
                m_hessian[N * (i1 + k) + j1 + l] = Hij_off;
                m_hessian[N * (i1 + l) + j1 + k] = Hij_off;
                m_hessian[N * (j1 + k) + i1 + l] = Hij_off;
                m_hessian[N * (j1 + l) + i1 + k] = Hij_off;
            }
        }
    }
};

/**
 * class which accumulates the energy one pair interaction at a time
 */
template <typename pairwise_interaction, typename distance_policy>
class NeighbourAccumulator {
    std::shared_ptr<pairwise_interaction> & m_interaction;
    std::shared_ptr<distance_policy> & m_dist;
    const size_t m_ndim;
    const double m_cutoff_sca;

public:
    pele::Array<std::vector<size_t>> m_neighbour_indss;
    pele::Array<std::vector<std::vector<double>>> m_neighbour_distss;

    NeighbourAccumulator(std::shared_ptr<pairwise_interaction> & interaction,
            std::shared_ptr<distance_policy> & dist,
            const size_t natoms, const size_t ndim, const double cutoff_sca)
        : m_interaction(interaction),
          m_dist(dist),
          m_ndim(ndim),
          m_cutoff_sca(cutoff_sca),
          m_neighbour_indss(natoms),
          m_neighbour_distss(natoms)
    {}

    void insert_atom_pair(Array<double> const & coords, const size_t atom_i, const size_t atom_j)
    {
        std::vector<double> dr(m_ndim);
        std::vector<double> neg_dr(m_ndim);
        m_dist->get_rij(dr.data(), coords.data() + m_ndim * atom_i, coords.data() + m_ndim * atom_j);
        double r2 = 0;
        for (size_t k = 0; k < m_ndim; ++k) {
            r2 += dr[k] * dr[k];
            neg_dr[k] = -dr[k];
        }
        const double r_H = m_interaction->m_radii[atom_i] + m_interaction->m_radii[atom_j];
        const double r_S = (1 + m_cutoff_sca) * r_H;
        const double r_S2 = r_S * r_S;
        if(r2 <= r_S2) {
            m_neighbour_indss[atom_i].push_back(atom_j);
            m_neighbour_indss[atom_j].push_back(atom_i);
            m_neighbour_distss[atom_i].push_back(dr);
            m_neighbour_distss[atom_j].push_back(neg_dr);
        }
    }
};

/**
 * Potential to loop over the list of atom pairs generated with the
 * cell list implementation in cell_lists.h.
 * This should also do the cell list construction and refresh, such that
 * the interface is the same for the user as with SimplePairwise.
 */
template <typename pairwise_interaction, typename distance_policy>
class CellListPotential : public PairwisePotentialInterface {
protected:
    const static size_t m_ndim = distance_policy::_ndim;
    pele::CellLists<distance_policy> m_cell_lists;
    std::shared_ptr<pairwise_interaction> m_interaction;
    std::shared_ptr<distance_policy> m_dist;
    const double m_radii_sca;
public:
    ~CellListPotential() {}
    CellListPotential(
            std::shared_ptr<pairwise_interaction> interaction,
            std::shared_ptr<distance_policy> dist,
            pele::Array<double> const & boxvec,
            double rcut, double ncellx_scale, const double radii_sca=0.0)
        : m_cell_lists(dist, boxvec, rcut, ncellx_scale),
          m_interaction(interaction),
          m_dist(dist),
          m_radii_sca(radii_sca)
    {}
    virtual size_t get_ndim(){return m_ndim;}

    virtual double get_energy(Array<double> & coords)
    {
        const size_t natoms = coords.size() / m_ndim;
        if (m_ndim * natoms != coords.size()) {
            throw std::runtime_error("coords.size() is not divisible by the number of dimensions");
        }

        update_iterator(coords);
        typedef EnergyAccumulator<pairwise_interaction, distance_policy> accumulator_t;
        accumulator_t accumulator(m_interaction, m_dist);
        auto looper = m_cell_lists.get_atom_pair_looper(accumulator);

        looper.loop_through_atom_pairs(coords);

        return accumulator.m_energy;
    }

    virtual double get_energy_gradient(Array<double> & coords, Array<double> & grad)
    {
        const size_t natoms = coords.size() / m_ndim;
        if (m_ndim * natoms != coords.size()) {
            throw std::runtime_error("coords.size() is not divisible by the number of dimensions");
        }
        if (coords.size() != grad.size()) {
            throw std::invalid_argument("the gradient has the wrong size");
        }

        update_iterator(coords);
        grad.assign(0.);
        typedef EnergyGradientAccumulator<pairwise_interaction, distance_policy> accumulator_t;
        accumulator_t accumulator(m_interaction, m_dist, grad);
        auto looper = m_cell_lists.get_atom_pair_looper(accumulator);

        looper.loop_through_atom_pairs(coords);

        return accumulator.m_energy;
    }

    virtual double get_energy_gradient_hessian(Array<double> & coords,
            Array<double> & grad, Array<double> & hess)
    {
        const size_t natoms = coords.size() / m_ndim;
        if (m_ndim * natoms != coords.size()) {
            throw std::runtime_error("coords.size() is not divisible by the number of dimensions");
        }
        if (coords.size() != grad.size()) {
            throw std::invalid_argument("the gradient has the wrong size");
        }
        if (hess.size() != coords.size() * coords.size()) {
            throw std::invalid_argument("the Hessian has the wrong size");
        }

        update_iterator(coords);
        grad.assign(0.);
        hess.assign(0.);
        typedef EnergyGradientHessianAccumulator<pairwise_interaction, distance_policy> accumulator_t;
        accumulator_t accumulator(m_interaction, m_dist, grad, hess);
        auto looper = m_cell_lists.get_atom_pair_looper(accumulator);

        looper.loop_through_atom_pairs(coords);

        return accumulator.m_energy;
    }

    virtual void get_neighbours(Array<double> & coords,
                                pele::Array<std::vector<size_t>> & neighbour_indss,
                                pele::Array<std::vector<std::vector<double>>> & neighbour_distss)
    {
        const size_t natoms = coords.size() / m_ndim;
        if (m_ndim * natoms != coords.size()) {
            throw std::runtime_error("coords.size() is not divisible by the number of dimensions");
        }
        if (m_interaction->m_radii.size() == 0) {
            throw std::runtime_error("Can't calculate neighbours, because the "
                                     "used interaction doesn't use radii. ");
        }

        update_iterator(coords);
        NeighbourAccumulator<pairwise_interaction, distance_policy> accumulator(
            m_interaction, m_dist, natoms, m_ndim, m_radii_sca);
        auto looper = m_cell_lists.get_atom_pair_looper(accumulator);

        looper.loop_through_atom_pairs(coords);

        neighbour_indss = accumulator.m_neighbour_indss;
        neighbour_distss = accumulator.m_neighbour_distss;
    }

    virtual inline size_t get_ndim() const { return m_ndim; }

    virtual inline void get_rij(double * const r_ij, double const * const r1, double const * const r2) const
    {
        return m_dist->get_rij(r_ij, r1, r2);
    }

    virtual inline double get_interaction_energy_gradient(double r2, double *gij, size_t atom_i, size_t atom_j) const
    {
        return m_interaction->energy_gradient(r2, gij, atom_i, atom_j);
    }

    virtual inline double get_interaction_energy_gradient_hessian(double r2, double *gij, double *hij, size_t atom_i, size_t atom_j) const
    {
        return m_interaction->energy_gradient_hessian(r2, gij, hij, atom_i, atom_j);
    }

protected:
    void update_iterator(Array<double> & coords)
    {
        m_cell_lists.update(coords);
    }
};

} //namespace pele

#endif //#ifndef _PELE_CELL_LIST_POTENTIAL_H
