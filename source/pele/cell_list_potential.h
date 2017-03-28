#ifndef _PELE_CELL_LIST_POTENTIAL_H
#define _PELE_CELL_LIST_POTENTIAL_H

#include <iostream>
#include <memory>
#include <vector>
#include <algorithm>
#include <utility>
#include <stdexcept>
#include <omp.h>

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
    std::vector<double*> m_energies;

public:
    ~EnergyAccumulator()
    {
        for(auto energy : m_energies) {
            delete energy;
        }
    }

    EnergyAccumulator(std::shared_ptr<pairwise_interaction> & interaction,
            std::shared_ptr<distance_policy> & dist)
        : m_interaction(interaction),
          m_dist(dist)
    {
        #ifdef _OPENMP
        m_energies = std::vector<double*>(omp_get_max_threads());
        #pragma omp parallel
        {
            m_energies[omp_get_thread_num()] = new double();
        }
        #else
        m_energies = std::vector<double*>(1);
        m_energies[0] = new double();
        #endif
    }

    void insert_atom_pair(Atom<m_ndim> const & atom_i, Atom<m_ndim> const & atom_j)
    {
        double dr[m_ndim];
        m_dist->get_rij(dr, atom_i.coords.data(), atom_j.coords.data());
        double r2 = 0;
        for (size_t k = 0; k < m_ndim; ++k) {
            r2 += dr[k] * dr[k];
        }
        const double radius_sum = atom_i.radius + atom_j.radius;
        #ifdef _OPENMP
        size_t tid = omp_get_thread_num();
        *m_energies[tid] += m_interaction->energy(r2, radius_sum);
        #else
        *m_energies[0] += m_interaction->energy(r2, radius_sum);
        #endif
    }

    double get_energy() {
        double energy = 0;
        for(size_t i = 0; i < m_energies.size(); i++) {
            energy += *m_energies[i];
        }
        return energy;
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
    std::vector<double*> m_energies;

public:
    ~EnergyGradientAccumulator()
    {
        for(auto energy : m_energies) {
            delete energy;
        }
    }

    EnergyGradientAccumulator(std::shared_ptr<pairwise_interaction> & interaction,
            std::shared_ptr<distance_policy> & dist)
        : m_interaction(interaction),
          m_dist(dist)
    {
        #ifdef _OPENMP
        m_energies = std::vector<double*>(omp_get_max_threads());
        #pragma omp parallel
        {
            m_energies[omp_get_thread_num()] = new double();
        }
        #else
        m_energies = std::vector<double*>(1);
        m_energies[0] = new double();
        #endif
    }

    void insert_atom_pair(Atom<m_ndim>& atom_i, Atom<m_ndim>& atom_j)
    {
        pele::VecN<m_ndim, double> dr;
        m_dist->get_rij(dr.data(), atom_i.coords.data(), atom_j.coords.data());
        double r2 = 0;
        for (size_t k = 0; k < m_ndim; ++k) {
            r2 += dr[k] * dr[k];
        }
        double gij;
        const double radius_sum = atom_i.radius + atom_j.radius;
        #ifdef _OPENMP
        size_t tid = omp_get_thread_num();
        *m_energies[tid] += m_interaction->energy_gradient(r2, &gij, radius_sum);
        #else
        *m_energies[0] += m_interaction->energy_gradient(r2, &gij, radius_sum);
        #endif
        if (gij != 0) {
            dr *= gij;
            atom_i.grad -= dr;
            atom_j.grad += dr;
        }
    }

    double get_energy() {
        double energy = 0;
        for(size_t i = 0; i < m_energies.size(); i++) {
            energy += *m_energies[i];
        }
        return energy;
    }

    void get_gradient(std::vector< std::vector<std::vector<pele::Atom<m_ndim>>>* > const & subdom_cell_atoms, pele::Array<double> & grad) {
        #ifdef _OPENMP
        #pragma omp parallel
        {
            size_t isubdom = omp_get_thread_num();
            get_gradient_subdom(*subdom_cell_atoms[isubdom], grad);
        }
        #else
        for (const auto cell_atoms : subdom_cell_atoms) {
            get_gradient_subdom(*cell_atoms, grad);
        }
        #endif
    }

    void get_gradient_subdom(std::vector<std::vector<pele::Atom<m_ndim>>> cell_atoms,
                             pele::Array<double> & grad) {
        for (auto const & atoms : cell_atoms) {
            for (auto const & atom : atoms) {
                const size_t ind_off = m_ndim * atom.index;
                for (size_t k = 0; k < m_ndim; k++) {
                    grad[ind_off + k] = atom.grad[k];
                }
            }
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
    std::vector<double*> m_energies;
    size_t m_natoms_ndim;

public:
    pele::Array<double> & m_hessian;

    ~EnergyGradientHessianAccumulator()
    {
        for(auto energy : m_energies) {
            delete energy;
        }
    }

    EnergyGradientHessianAccumulator(std::shared_ptr<pairwise_interaction> & interaction,
            std::shared_ptr<distance_policy> & dist,
            pele::Array<double> & hessian, const size_t natoms_ndim)
        : m_interaction(interaction),
          m_dist(dist),
          m_hessian(hessian),
          m_natoms_ndim(natoms_ndim)
    {
        #ifdef _OPENMP
        m_energies = std::vector<double*>(omp_get_max_threads());
        #pragma omp parallel
        {
            m_energies[omp_get_thread_num()] = new double();
        }
        #else
        m_energies = std::vector<double*>(1);
        m_energies[0] = new double();
        #endif
    }

    void insert_atom_pair(Atom<m_ndim> & atom_i, Atom<m_ndim> & atom_j)
    {
        pele::VecN<m_ndim, double> dr;
        m_dist->get_rij(dr.data(), atom_i.coords.data(), atom_j.coords.data());
        double r2 = 0;
        for (size_t k = 0; k < m_ndim; ++k) {
            r2 += dr[k] * dr[k];
        }
        double gij, hij;
        const double radius_sum = atom_i.radius + atom_j.radius;
        #ifdef _OPENMP
        size_t tid = omp_get_thread_num();
        *m_energies[tid] += m_interaction->energy_gradient_hessian(r2, &gij, &hij, radius_sum);
        #else
        *m_energies[0] += m_interaction->energy_gradient_hessian(r2, &gij, &hij, radius_sum);
        #endif
        if (gij != 0) {
            atom_i.grad -= dr * gij;
            atom_j.grad += dr * gij;
        }
        //this part is copied from simple_pairwise_potential.h
        //(even more so than the rest)
        const size_t N = m_natoms_ndim;
        const size_t i1 = m_ndim * atom_i.index;
        const size_t j1 = m_ndim * atom_j.index;
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

    double get_energy() {
        double energy = 0;
        for(size_t i = 0; i < m_energies.size(); i++) {
            energy += *m_energies[i];
        }
        return energy;
    }

    void get_gradient(std::vector< std::vector<std::vector<pele::Atom<m_ndim>>>* > const & subdom_cell_atoms, pele::Array<double> & grad) {
        #ifdef _OPENMP
        #pragma omp parallel
        {
            size_t isubdom = omp_get_thread_num();
            get_gradient_subdom(*subdom_cell_atoms[isubdom], grad);
        }
        #else
        for (const auto cell_atoms : subdom_cell_atoms) {
            get_gradient_subdom(*cell_atoms, grad);
        }
        #endif
    }

    void get_gradient_subdom(std::vector<std::vector<pele::Atom<m_ndim>>> cell_atoms,
                             pele::Array<double> & grad) {
        for (auto const & atoms : cell_atoms) {
            for (auto const & atom : atoms) {
                const size_t ind_off = m_ndim * atom.index;
                for (size_t k = 0; k < m_ndim; k++) {
                    grad[ind_off + k] = atom.grad[k];
                }
            }
        }
    }
};

/**
 * class which accumulates the energy one pair interaction at a time
 */
template <typename pairwise_interaction, typename distance_policy>
class NeighborAccumulator {
    std::shared_ptr<pairwise_interaction> & m_interaction;
    std::shared_ptr<distance_policy> & m_dist;
    const static size_t m_ndim = distance_policy::_ndim;
    const double m_cutoff_sca;
    pele::Array<short> const & m_include_atoms;

public:
    pele::Array<std::vector<size_t>> m_neighbor_indss;
    pele::Array<std::vector<std::vector<double>>> m_neighbor_distss;

    NeighborAccumulator(std::shared_ptr<pairwise_interaction> & interaction,
            std::shared_ptr<distance_policy> & dist,
            const size_t natoms, const double cutoff_sca,
            pele::Array<short> const & include_atoms)
        : m_interaction(interaction),
          m_dist(dist),
          m_cutoff_sca(cutoff_sca),
          m_include_atoms(include_atoms),
          m_neighbor_indss(natoms),
          m_neighbor_distss(natoms)
    {}

    void insert_atom_pair(Atom<m_ndim> const & atom_i, Atom<m_ndim> const & atom_j)
    {
        if (m_include_atoms[atom_i.index] && m_include_atoms[atom_j.index]) {
            std::vector<double> dr(m_ndim);
            std::vector<double> neg_dr(m_ndim);
            m_dist->get_rij(dr.data(), atom_i.coords.data(), atom_j.coords.data());
            double r2 = 0;
            for (size_t k = 0; k < m_ndim; ++k) {
                r2 += dr[k] * dr[k];
                neg_dr[k] = -dr[k];
            }
            const double radius_sum = atom_i.radius + atom_j.radius;
            const double r_S = m_cutoff_sca * radius_sum;
            const double r_S2 = r_S * r_S;
            if(r2 <= r_S2) {
                m_neighbor_indss[atom_i.index].push_back(atom_j.index);
                m_neighbor_indss[atom_j.index].push_back(atom_i.index);
                m_neighbor_distss[atom_i.index].push_back(dr);
                m_neighbor_distss[atom_j.index].push_back(neg_dr);
            }
        }
    }
};

/**
 * class which accumulates the energy one pair interaction at a time
 */
template <typename pairwise_interaction, typename distance_policy>
class OverlapAccumulator {
    std::shared_ptr<pairwise_interaction> & m_interaction;
    std::shared_ptr<distance_policy> & m_dist;
    const static size_t m_ndim = distance_policy::_ndim;

public:
    std::vector<size_t> m_overlap_inds;

    OverlapAccumulator(std::shared_ptr<pairwise_interaction> & interaction,
            std::shared_ptr<distance_policy> & dist,
            const size_t natoms)
        : m_interaction(interaction),
          m_dist(dist)
    {}

    void insert_atom_pair(Atom<m_ndim> const & atom_i, Atom<m_ndim> const & atom_j)
    {
        std::vector<double> dr(m_ndim);
        m_dist->get_rij(dr.data(), atom_i.coords.data(), atom_j.coords.data());
        double r2 = 0;
        for (size_t k = 0; k < m_ndim; ++k) {
            r2 += dr[k] * dr[k];
        }
        const double radius_sum = atom_i.radius + atom_j.radius;
        const double r_H2 = radius_sum * radius_sum;
        if(r2 <= r_H2) {
            #ifdef _OPENMP
            #pragma omp critical
            {
                m_overlap_inds.push_back(atom_i.index);
                m_overlap_inds.push_back(atom_j.index);
            }
            #else
            m_overlap_inds.push_back(atom_i.index);
            m_overlap_inds.push_back(atom_j.index);
            #endif
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
            double rcut, double ncellx_scale,
            const pele::Array<double> radii,
            const double radii_sca=0.0)
        : PairwisePotentialInterface(radii),
          m_cell_lists(dist, boxvec, rcut, ncellx_scale, radii),
          m_interaction(interaction),
          m_dist(dist),
          m_radii_sca(radii_sca)
    {}

    CellListPotential(
            std::shared_ptr<pairwise_interaction> interaction,
            std::shared_ptr<distance_policy> dist,
            pele::Array<double> const & boxvec,
            double rcut, double ncellx_scale)
        : m_cell_lists(dist, boxvec, rcut, ncellx_scale),
          m_interaction(interaction),
          m_dist(dist),
          m_radii_sca(0.0)
    {}

    virtual size_t get_ndim(){return m_ndim;}

    virtual double get_energy(Array<double> const & coords)
    {
        const size_t natoms = coords.size() / m_ndim;
        if (m_ndim * natoms != coords.size()) {
            throw std::runtime_error("coords.size() is not divisible by the number of dimensions");
        }

        update_iterator(coords);
        typedef EnergyAccumulator<pairwise_interaction, distance_policy> accumulator_t;
        accumulator_t accumulator(m_interaction, m_dist);
        auto looper = m_cell_lists.get_atom_pair_looper(accumulator);

        looper.loop_through_atom_pairs();

        return accumulator.get_energy();
    }

    virtual double get_energy_gradient(Array<double> const & coords, Array<double> & grad)
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
        accumulator_t accumulator(m_interaction, m_dist);
        auto looper = m_cell_lists.get_atom_pair_looper(accumulator);

        looper.loop_through_atom_pairs();

        accumulator.get_gradient(m_cell_lists.get_atoms(), grad);
        return accumulator.get_energy();
    }

    virtual double get_energy_gradient_hessian(Array<double> const & coords,
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
        accumulator_t accumulator(m_interaction, m_dist, hess, coords.size());
        auto looper = m_cell_lists.get_atom_pair_looper(accumulator);

        looper.loop_through_atom_pairs();

        accumulator.get_gradient(m_cell_lists.get_atoms(), grad);
        return accumulator.get_energy();
    }

    virtual void get_neighbors(pele::Array<double> const & coords,
                                pele::Array<std::vector<size_t>> & neighbor_indss,
                                pele::Array<std::vector<std::vector<double>>> & neighbor_distss,
                                const double cutoff_factor = 1.0)
    {
        size_t natoms = coords.size() / m_ndim;
        pele::Array<short> include_atoms(natoms, 1);
        get_neighbors_picky(coords, neighbor_indss, neighbor_distss, include_atoms, cutoff_factor);
    }

    virtual void get_neighbors_picky(pele::Array<double> const & coords,
                                      pele::Array<std::vector<size_t>> & neighbor_indss,
                                      pele::Array<std::vector<std::vector<double>>> & neighbor_distss,
                                      pele::Array<short> const & include_atoms,
                                      const double cutoff_factor = 1.0)
    {
        const size_t natoms = coords.size() / m_ndim;
        if (m_ndim * natoms != coords.size()) {
            throw std::runtime_error("coords.size() is not divisible by the number of dimensions");
        }
        if (natoms != include_atoms.size()) {
            throw std::runtime_error("include_atoms.size() is not equal to the number of atoms");
        }
        if (m_radii.size() == 0) {
            throw std::runtime_error("Can't calculate neighbors, because the "
                                     "used interaction doesn't use radii. ");
        }

        update_iterator(coords);
        NeighborAccumulator<pairwise_interaction, distance_policy> accumulator(
            m_interaction, m_dist, natoms, (1 + m_radii_sca) * cutoff_factor, include_atoms);
        auto looper = m_cell_lists.get_atom_pair_looper(accumulator);

        looper.loop_through_atom_pairs();

        neighbor_indss = accumulator.m_neighbor_indss;
        neighbor_distss = accumulator.m_neighbor_distss;
    }

    virtual std::vector<size_t> get_overlaps(Array<double> const & coords)
    {
        const size_t natoms = coords.size() / m_ndim;
        if (m_ndim * natoms != coords.size()) {
            throw std::runtime_error("coords.size() is not divisible by the number of dimensions");
        }
        if (m_radii.size() == 0) {
            throw std::runtime_error("Can't calculate neighbors, because the "
                                     "used interaction doesn't use radii. ");
        }

        update_iterator(coords);
        OverlapAccumulator<pairwise_interaction, distance_policy> accumulator(
            m_interaction, m_dist, natoms);
        auto looper = m_cell_lists.get_atom_pair_looper(accumulator);

        looper.loop_through_atom_pairs();

        return accumulator.m_overlap_inds;
    }

    virtual inline size_t get_ndim() const { return m_ndim; }

    virtual inline void get_rij(double * const r_ij, double const * const r1, double const * const r2) const
    {
        return m_dist->get_rij(r_ij, r1, r2);
    }

    virtual inline double get_interaction_energy_gradient(double r2, double *gij, size_t atom_i, size_t atom_j) const
    {
        return m_interaction->energy_gradient(r2, gij, sum_radii(atom_i, atom_j));
    }

    virtual inline double get_interaction_energy_gradient_hessian(double r2, double *gij, double *hij, size_t atom_i, size_t atom_j) const
    {
        return m_interaction->energy_gradient_hessian(r2, gij, hij, sum_radii(atom_i, atom_j));
    }

protected:
    void update_iterator(Array<double> const & coords)
    {
        m_cell_lists.update(coords);
    }
};

} //namespace pele

#endif //#ifndef _PELE_CELL_LIST_POTENTIAL_H
