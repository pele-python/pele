#ifndef _PELE_CELL_LISTS_H_
#define _PELE_CELL_LISTS_H_

#include <iostream>
#include <memory>
#include <exception>
#include <cassert>
#include <vector>
#include <algorithm>
#include <omp.h>
#include <sstream>
#include <array>
#include <limits>

#include "base_potential.h"
#include "array.h"
#include "distance.h"
#include "vecn.h"
#include "queue.h"

namespace {

static const long CELL_END = -1;

template<class T, size_t box_dimension>
struct periodic_policy_check_helper {
    const static bool is_periodic = false;
};

template<size_t box_dimension>
struct periodic_policy_check_helper<pele::periodic_distance<box_dimension>, box_dimension > {
    const static bool is_periodic = true;
};

template<class T>
struct periodic_policy_check {
    const static bool is_periodic = periodic_policy_check_helper<T, T::_ndim>::is_periodic;
};

/**
 * this iterator facilitates looping through the atoms in a cell.
 */
template<size_t ndim>
class AtomInCellIterator {
private:
    std::vector<long> const * m_atoms;
    long m_current_atom;

public:
    AtomInCellIterator(std::vector<long> const * atoms, long first_atom)
        : m_atoms(atoms),
          m_current_atom(first_atom)
    {}

    /**
     * return the current atom
     */
    inline long operator*() const
    {
        return m_current_atom;
    }

    /**
     * advance to the next atom
     *
     * This will seg fault if we go past the end of the list.  It
     * is up to the user to ensure that doesn't happen.
     */
    inline void operator++()
    {
        m_current_atom = (*m_atoms)[m_current_atom];
    }
};

/**
 * container for the cell lists
 */
template<size_t ndim>
class CellListsContainer {
protected:
    /** Construct m_cells for each subdomain
     */
    void setup_cells(const size_t nsubdoms, const std::vector<size_t> subdom_ncells) {
        #ifdef _OPENMP
        #pragma omp parallel
    	{
            size_t isubdom = omp_get_thread_num();
            m_cells[isubdom] = std::vector<long>(subdom_ncells[isubdom], CELL_END);
    	}
        #else
        for (size_t isubdom = 0; isubdom < nsubdoms; ++isubdom) {
            m_cells[isubdom] = std::vector<long>(subdom_ncells[isubdom], CELL_END);
    	}
        #endif
    }

public:
    /**
     * m_atoms is a vector of atom indices.
     *
     * Each atom index points to the next atom index in the same cell.
     * If there are no more atoms in this cell, the value is set to CELL_END.
     */
    std::vector<long> m_atoms;

    /**
     * m_cells is a vector of vectors of atom indices.
     *
     * Each atom index corresponds to the first atom in this cell. The atom index also
     * serves as index in m_atoms, where the next atom index in this cell is saved.
     * If there are no atoms in this cell, the value is set to CELL_END.
     */
    std::vector< std::vector<long> > m_cells;

    /** vectors of pairs of neighboring cells inside each subdomain*/
    std::vector< std::vector< std::array<long*, 2> > > m_cell_neighbor_pairs_inner;

    /** vectors of pairs of neighboring cells between subdomains*/
    std::vector< std::vector< std::array<long*, 2> > > m_cell_neighbor_pairs_boundary;

    CellListsContainer(const std::vector<size_t> subdom_ncells)
        : m_cells(subdom_ncells.size()),
          m_cell_neighbor_pairs_inner(subdom_ncells.size()),
          m_cell_neighbor_pairs_boundary(subdom_ncells.size())
    {
        setup_cells(subdom_ncells.size(), subdom_ncells);
    }

    /**
     * reset all data in the vectors
     */
    void reset()
    {
        m_atoms.assign(m_atoms.size(), CELL_END);
        for (auto subdom_cells : m_cells) {
            subdom_cells.assign(subdom_cells.size(), CELL_END);
        }
    }

    /**
     * add an atom to a cell
     */
    inline void add_atom_to_cell(const long iatom, const size_t icell, const size_t isubdom)
    {
        m_atoms[iatom] = m_cells[isubdom][icell];
        m_cells[isubdom][icell] = iatom;
    }

    /**
     * remove an atom from a cell using its and its parent's index
     */
    inline void remove_atom_from_cell(const long parent_index, const long index,
                                      const size_t icell, const size_t isubdom)
    {
        if (parent_index == CELL_END) {
            // Remove first atom in cell
            m_cells[isubdom][icell] = m_atoms[index];
        } else {
            // Remove any other atom in cell
            m_atoms[parent_index] = m_atoms[index];
        }
    }

    /**
     * remove an atom from a cell using only its own index
     */
    inline void remove_atom_from_cell(const long index, const size_t icell, const size_t isubdom)
    {
        long parent_index = CELL_END;
        for (auto cell_iter = getIterator(icell, isubdom);
             *cell_iter != CELL_END;
             ++cell_iter) {
            if (*cell_iter == index) {
                remove_atom_from_cell(parent_index, index, icell, isubdom);
                break;
            } else {
                parent_index = *cell_iter;
            }
        }
    }

    /**
     * return an iterator over the atoms in cell `icell` of subdomain `isubdom`
     */
    AtomInCellIterator<ndim> getIterator(const size_t icell, const size_t isubdom) const
    {
        return AtomInCellIterator<ndim>(&m_atoms, m_cells[isubdom][icell]);
    }

    /**
     * return an iterator over the atoms in cell
     */
    AtomInCellIterator<ndim> getIterator(const long* cell_ptr) const
    {
        return AtomInCellIterator<ndim>(&m_atoms, *cell_ptr);
    }
};

}

namespace pele {

/**
 * this does the looping over atom pairs within the cell lists framework.
 *
 * This uses the visitor design pattern, in the same way that the Boost graph
 * uses visitors for, e.g., breadth_first_search.
 * http://www.boost.org/doc/libs/1_56_0/boost/graph/breadth_first_search.hpp
 *
 * The visitor is called for every pair of atoms in the system.  It
 * is meant to be used for, e.g. accumulating the energy or the gradient.
 *
 * We use a template rather than an interface because the visitor will be called many
 * times over a short period and the additional overhead of an interface might be a problem.
 */
template <class visitor_t, size_t ndim>
class CellListsLoop {
protected:
    visitor_t & m_visitor;
    CellListsContainer<ndim> const & m_container;

    void loop_cell_pairs(
        std::vector< std::array<long*, 2> >
        const & neighbor_pairs, size_t isubdom)
    {
        for (auto const & ijpair : neighbor_pairs) {
            // do double loop through atoms, avoiding duplicate pairs
            for (auto icell_iter = m_container.getIterator(ijpair[0]);
                 *icell_iter != CELL_END;
                 ++icell_iter) {
                // if icell==jcell we need to avoid duplicate atom pairs
                auto jend = (ijpair[0] == ijpair[1]) ? *icell_iter : CELL_END;
                for (auto jcell_iter = m_container.getIterator(ijpair[1]);
                     *jcell_iter != jend;
                     ++jcell_iter) {
                    m_visitor.insert_atom_pair(*icell_iter, *jcell_iter, isubdom);
                }
            }
        }
    }

public:
    CellListsLoop(visitor_t & visitor, CellListsContainer<ndim> const & container)
        : m_visitor(visitor),
          m_container(container)
    {}

    void loop_through_atom_pairs()
    {
        #ifdef _OPENMP
        #pragma omp parallel
        {
            size_t isubdom = omp_get_thread_num();
            loop_cell_pairs(m_container.m_cell_neighbor_pairs_inner[isubdom], isubdom);
            #pragma omp barrier
            loop_cell_pairs(m_container.m_cell_neighbor_pairs_boundary[isubdom], isubdom);
        }
        #else
        size_t nsubdoms = m_container.m_cells.size();
        for(size_t isubdom = 0; isubdom < nsubdoms; ++isubdom) {
            loop_cell_pairs(m_container.m_cell_neighbor_pairs_inner[isubdom], isubdom);
            loop_cell_pairs(m_container.m_cell_neighbor_pairs_boundary[isubdom], isubdom);
        }
        #endif
    }
};

template<typename distance_policy>
class LatticeNeighbors {
protected:
    /** Calculate number of cells of each subdomain in the split direction (y)
     */
    void calc_subdom_stats() {
        size_t remainder = m_ncells_vec[1] % m_nsubdoms;
        m_subdoms_balanced = remainder == 0;
        size_t subdom_len;
        m_subdom_avg_len = 0;
        m_subdom_ncells = std::vector<size_t>(m_nsubdoms);
        m_subdom_limits = std::vector<size_t>(m_nsubdoms + 1);
        m_subdom_limits[0] = 0;
        std::stringstream ncells_stream;
        ncells_stream << "Length per subdomain: ";
        for (size_t isubdom = 0; isubdom < m_nsubdoms; ++isubdom) {
            if(remainder > isubdom) {
                subdom_len = m_ncells_vec[1] / m_nsubdoms + 1;
            } else {
                subdom_len = m_ncells_vec[1] / m_nsubdoms;
            }
            m_subdom_avg_len += subdom_len;
            m_subdom_limits[isubdom + 1] = m_subdom_limits[isubdom] + subdom_len;
            m_subdom_ncells[isubdom] = subdom_len * m_ncells / m_ncells_vec[1];
            ncells_stream << subdom_len;
            if(isubdom < m_nsubdoms - 1) {
                 ncells_stream << ", ";
            }
        }
        m_subdom_avg_len /= m_nsubdoms;
        std::cout << ncells_stream.str() << std::endl;
        #ifdef _OPENMP
        if( (m_boxvec[1] / m_rcut) / m_nsubdoms < 2 && omp_get_max_threads() > 1) {
            throw std::runtime_error("Some subdomains are thinner than 2*rcut. The "
                "parallelization would not be stable. Reduce the number of threads "
                "(environment variable OMP_NUM_THREADS)!");
        }
        #else
        if( m_ncells_vec[1] / m_nsubdoms < 1) {
            throw std::runtime_error("Some subdomains are thinner than 1 row!");
        }
        #endif
    }

public:
    static const size_t ndim = distance_policy::_ndim;
    const std::shared_ptr<distance_policy> m_dist; //!< the distance function

    typedef VecN<ndim, size_t> cell_vec_t;
    pele::VecN<ndim> m_boxvec;
    pele::VecN<ndim> m_inv_boxvec; //!< inverse of boxvec
    double m_rcut;
    cell_vec_t m_ncells_vec; //!< the number of cells in each dimension
    pele::VecN<ndim> m_rcell_vec; //!< the cell length in each dimension
    size_t m_ncells;
    size_t m_nsubdoms;
    std::vector<size_t> m_subdom_limits; //!< boundary indices of each subdomain in the split direction (y)
    std::vector<size_t> m_subdom_ncells; //!< number of cells of each subdomain
    bool m_subdoms_balanced; //!< true if all subdomains have the same number of cells
    size_t m_subdom_avg_len; //!< Average subdomain length. Only used in balanced case, therefore integer.

    LatticeNeighbors(std::shared_ptr<distance_policy> const & dist,
            pele::Array<double> const & boxvec,
            const double rcut,
            pele::Array<size_t> const & ncells_vec)
        : m_dist(dist),
          m_boxvec(boxvec),
          m_rcut(rcut),
          m_ncells_vec(ncells_vec),
          m_inv_boxvec(ndim)
    {
        #ifdef _OPENMP
        m_nsubdoms = omp_get_max_threads();
        #else
        m_nsubdoms = 1;
        #endif

        m_ncells = m_ncells_vec.prod();
        for (size_t idim = 0; idim < ndim; ++idim) {
            m_inv_boxvec[idim] = 1 / m_boxvec[idim];
            m_rcell_vec[idim] = m_boxvec[idim] / m_ncells_vec[idim];
        }
        calc_subdom_stats();
    }

    /**
     * convert a cell vector to the global cell index
     * The cell vector needs to be within range [0, m_ncells_vec]
     */
    size_t cell_vec_to_global_ind(cell_vec_t const & v) const
    {
        // note: the speed of this function is important because it is called
        // once per atom each time get_energy is called.
        size_t cum = 1;
        size_t index = 0;
        for (size_t idim = 0; idim < ndim; ++idim) {
            index += v[idim] * cum;
            cum *= m_ncells_vec[idim];
        }
        return index;
    }

    /**
     * convert a cell vector to the local cell index
     * The cell vector needs to be within range [0, m_ncells_vec]
     */
    size_t cell_vec_to_local_ind(cell_vec_t const & v, const size_t isubdom) const
    {
        // note: the speed of this function is important because it is called
        // once per atom each time get_energy is called.
        size_t cum = 1;
        size_t index = 0;
        for (size_t idim = 0; idim < ndim; ++idim) {
            if (idim == 1) {
                index += (v[1] - m_subdom_limits[isubdom]) * cum;
                cum *= m_subdom_limits[isubdom + 1] - m_subdom_limits[isubdom];
            } else {
                index += v[idim] * cum;
                cum *= m_ncells_vec[idim];
            }
        }
        return index;
    }

    /**
     * convert a global cell index to a cell vector
     */
    cell_vec_t global_ind_to_cell_vec(const size_t icell) const
    {
        cell_vec_t v;
        size_t remaining_icell = icell;
        size_t fraction;
        for (size_t idim = 0; idim < ndim - 1; ++idim) {
            fraction = (size_t) remaining_icell / m_ncells_vec[idim];
            // should be optimized by compiler to use remainder of previous division
            v[idim] = remaining_icell % m_ncells_vec[idim];
            remaining_icell = fraction;
        }
        v[ndim - 1] = remaining_icell;
        return v;
    }

    /**
     * convert a local cell index to a cell vector
     */
    cell_vec_t local_ind_to_cell_vec(const size_t icell, const size_t isubdom) const
    {
        cell_vec_t v;
        size_t remaining_icell = icell;
        size_t fraction;
        for (size_t idim = 0; idim < ndim - 1; ++idim) {
            if (idim == 1) {
                size_t split_len = m_subdom_limits[isubdom + 1] - m_subdom_limits[isubdom];
                fraction = (size_t) remaining_icell / split_len;
                // should be optimized by compiler to use remainder of previous division
                v[idim] = (remaining_icell % split_len) + m_subdom_limits[isubdom];
            } else {
                fraction = (size_t) remaining_icell / m_ncells_vec[idim];
                // should be optimized by compiler to use remainder of previous division
                v[idim] = remaining_icell % m_ncells_vec[idim];
            }
            remaining_icell = fraction;
        }
        v[ndim - 1] = remaining_icell;
        if (ndim - 1 == 1) {
            v[ndim - 1] += m_subdom_limits[isubdom];
        }
        return v;
    }

    /**
     * convert a cell vector to a positional vector in real space
     */
    VecN<ndim> cell_vec_to_position(cell_vec_t const & v) const
    {
        VecN<ndim> x;
        for (size_t idim = 0; idim < ndim; ++idim) {
            x[idim] = m_rcell_vec[idim] * v[idim];
        }
        return x;
    }

    /**
     * return the cell vector of the cell containing position x
     */
    cell_vec_t position_to_cell_vec(const double * const x) const
    {
        double x_in_box[ndim];
        std::copy(x, x + ndim, x_in_box);
        m_dist->put_atom_in_box(x_in_box);

        cell_vec_t cell_vec;
        for(size_t idim = 0; idim < ndim; ++idim) {
            cell_vec[idim] = m_ncells_vec[idim] * (x_in_box[idim] * m_inv_boxvec[idim] + 0.5
                                                   - std::numeric_limits<double>::epsilon());
            assert(cell_vec[idim] < m_ncells_vec[idim]); // Cell index is inside bounds (if not, the coordinates might be defective/too large)
        }
        return cell_vec;
    }

    /**
     * return the subdomain and the local index of the cell that contains position x
     */
    void position_to_local_ind(const double * const x, size_t & icell, size_t & isubdom) const
    {
        // note: the speed of this function is important because it is called
        // once per atom each time update_container is called.
        cell_vec_t cell_vec = position_to_cell_vec(x);
        isubdom = get_subdomain(cell_vec);
        icell = cell_vec_to_local_ind(cell_vec, isubdom);
    }

    /**
     * convert a global cell index to the subdomain and local index
     */
    void global_ind_to_local_ind(const size_t icell_global, size_t & icell_local, size_t & isubdom) const
    {
        cell_vec_t cell_vec = global_ind_to_cell_vec(icell_global);
        isubdom = get_subdomain(cell_vec);
        icell_local = cell_vec_to_local_ind(cell_vec, isubdom);
    }

    /**
     * convert local cell index and subdomain to the global index
     */
    size_t local_ind_to_global_ind(const size_t icell_local, const size_t isubdom) const
    {
        cell_vec_t cell_vec = local_ind_to_cell_vec(icell_local, isubdom);
        return cell_vec_to_global_ind(cell_vec);
    }

    /**
     * return the subdomain containing the cell
     */
    size_t get_subdomain(cell_vec_t const & v) const
    {
        if(m_subdoms_balanced) {
            return v[1] / m_subdom_avg_len;
        } else {
            return std::upper_bound(m_subdom_limits.begin(), m_subdom_limits.end(), v[1])
                   - m_subdom_limits.begin() - 1;
        }
    }

    /**
     * return the minimum corner to corner distance between two cells
     */
    double minimum_distance(cell_vec_t const & v1, cell_vec_t const & v2) const
    {
        // copy them so we don't accidentally change them
        auto lower_left1 = cell_vec_to_position(v1);
        auto lower_left2 = cell_vec_to_position(v2);
        pele::VecN<ndim> ll1, ll2, dr;
        pele::VecN<ndim> minimum_dist; // the minimum possible distance in each direction
        for (size_t i = 0; i < ndim; ++i) {
            double min_dist = std::numeric_limits<double>::max();
            double dri;
            double rcell = m_rcell_vec[i];
            // find the minimum distance in the i'th direction.
            ll1 = lower_left1;
            ll2 = lower_left2;
            m_dist->get_rij(dr.data(), ll1.data(), ll2.data());
            dri = std::abs(dr[i]);
            if (dri < min_dist) {
                min_dist = dri;
            }

            ll1 = lower_left1;
            ll2 = lower_left2;
            ll1[i] += rcell;
            m_dist->get_rij(dr.data(), ll1.data(), ll2.data());
            dri = std::abs(dr[i]);
            if (dri < min_dist) {
                min_dist = dri;
            }

            ll1 = lower_left1;
            ll2 = lower_left2;
            ll2[i] += rcell;
            m_dist->get_rij(dr.data(), ll1.data(), ll2.data());
            dri = std::abs(dr[i]);
            if (dri < min_dist) {
                min_dist = dri;
            }

            ll1 = lower_left1;
            ll2 = lower_left2;
            ll1[i] += rcell;
            ll2[i] += rcell;
            m_dist->get_rij(dr.data(), ll1.data(), ll2.data());
            dri = std::abs(dr[i]);
            if (dri < min_dist) {
                min_dist = dri;
            }

            minimum_dist[i] = min_dist;
        }
        double r2_min = dot<ndim> (minimum_dist, minimum_dist);
        return std::sqrt(r2_min);
    }

    /**
     * recursive function to find the neighbors of a given cell
     * using global cell indices
     */
    void find_global_neighbor_inds(const size_t idim, cell_vec_t const & v0,
            std::vector<size_t> & neighbors,
            cell_vec_t const & vorigin
            ) const
    {
        if (idim == ndim) {
            double rmin = minimum_distance(v0, vorigin);
            if (rmin <= m_rcut) {
                neighbors.push_back(cell_vec_to_global_ind(v0));
            }
        } else {
            auto v = v0;
            for(v[idim] = 0; v[idim] < m_ncells_vec[idim]; ++v[idim]) {
                find_global_neighbor_inds(idim+1, v, neighbors, vorigin);
            }
        }
    }

    /**
     * return a vector of all the neighbors of icell (including icell itself)
     * using global cell indices
     */
    std::vector<size_t> find_all_global_neighbor_inds(const size_t icell) const
    {
        auto vcell = global_ind_to_cell_vec(icell);

        std::vector<size_t> neighbors;
        find_global_neighbor_inds(0, vcell, neighbors, vcell);
        return neighbors;
    }

    /**
     * return a list of all pairs of neighboring cells
     */
    void find_neighbor_pairs(
        std::vector< std::vector< std::array<long*, 2> > > & cell_neighbors_inner,
        std::vector< std::vector< std::array<long*, 2> > > & cell_neighbors_boundary,
        std::vector< std::vector<long> > & cells) const
    {
        #ifdef _OPENMP
        #pragma omp parallel
        {
            size_t isubdom = omp_get_thread_num();
            find_neighbor_pairs_subdom(
                cell_neighbors_inner, cell_neighbors_boundary, cells, isubdom);
        }
        #else
        for(size_t isubdom = 0; isubdom < m_nsubdoms; ++isubdom) {
            find_neighbor_pairs_subdom(
                cell_neighbors_inner, cell_neighbors_boundary, cells, isubdom);
        }
        #endif
    }

    /**
     * return a list of all pairs of neighboring cells originating from a subdomain
     */
    void find_neighbor_pairs_subdom(
        std::vector< std::vector< std::array<long*, 2> > > & cell_neighbors_inner,
        std::vector< std::vector< std::array<long*, 2> > > & cell_neighbors_boundary,
        std::vector< std::vector<long> > & cells,
        const size_t isubdom) const
    {
        cell_neighbors_inner[isubdom] = std::vector< std::array<long*, 2> >();
        cell_neighbors_boundary[isubdom] = std::vector< std::array<long*, 2> >();

        // Reserve memory for the cell neighbors (only for performance)
        // This calculates the exact number of cells if ncellx_scale <= 1
        if(m_nsubdoms == 1) {
            cell_neighbors_inner[isubdom].reserve(m_subdom_ncells[isubdom] * ((std::pow(3, ndim) + 1) * 0.5));
        } else {
            const size_t subdom_length = m_subdom_limits[isubdom + 1] - m_subdom_limits[isubdom];
            const size_t surface_ncells = m_subdom_ncells[isubdom] / subdom_length;
            if(isubdom == 0) {
                // For Lees-Edwards Boundary Conditions
                cell_neighbors_boundary[isubdom].reserve(surface_ncells * 4 * std::pow(3, ndim - 2));
            } else {
                cell_neighbors_boundary[isubdom].reserve(surface_ncells * std::pow(3, ndim - 1));
            }
            cell_neighbors_inner[isubdom].reserve(m_subdom_ncells[isubdom] * ((std::pow(3, ndim) + 1) * 0.5) - surface_ncells * std::pow(3, ndim - 1));
        }

        for (size_t local_icell = 0; local_icell < m_subdom_ncells[isubdom]; ++local_icell) {
            size_t global_icell = local_ind_to_global_ind(local_icell, isubdom);
            auto global_neighbors = find_all_global_neighbor_inds(global_icell);
            for (size_t global_jcell : global_neighbors) {
                size_t local_jcell, jsubdom;
                global_ind_to_local_ind(global_jcell, local_jcell, jsubdom);
                if(isubdom == jsubdom)
                {
                    if (local_jcell >= local_icell) { // avoid duplicates
                        std::array<long*, 2> neighbors = {&cells[isubdom][local_icell],
                                                            &cells[jsubdom][local_jcell]};
                        cell_neighbors_inner[isubdom].push_back(neighbors);
                    }
                } else {
                    if(pos_direction_y(global_icell, global_jcell)) { // avoid duplicates, balance load
                        std::array<long*, 2> neighbors = {&cells[isubdom][local_icell],
                                                            &cells[jsubdom][local_jcell]};
                        cell_neighbors_boundary[isubdom].push_back(neighbors);
                    }
                }
            }
        }
    }

    /** Check if the direction of this neighborhood is positive
     *
     * True, if the direction is downward facing considering boundary conditions
     */
    bool pos_direction_y(const size_t icell, const size_t jcell) const {
        auto ipos = cell_vec_to_position(global_ind_to_cell_vec(icell));
        auto jpos = cell_vec_to_position(global_ind_to_cell_vec(jcell));
        pele::VecN<ndim> dr;
        m_dist->get_rij(dr.data(), ipos.data(), jpos.data());
        return dr[1] >= 0;
    }
};

/**
 * cell list currently only work with box of equal side lengths
 * cell lists are currently not implemented for non cubic boxes:
 * in that case m_ncellx should be an array rather than a scalar and the definition
 * of ncells and rcell would change to array. This implies that in all the function these
 * would need to be replace with the correct array element. This adds room for errors
 * so in this first implementation we do not account for that scenario
 */
template<typename distance_policy = periodic_distance<3> >
class CellLists{
public:
    static const size_t m_ndim = distance_policy::_ndim;
protected:

    bool m_initialized; // flag for whether the cell lists have been initialized with coordinates
    pele::LatticeNeighbors<distance_policy> m_lattice_tool;
    std::vector<SafePushQueue<std::array<long, 2>>> add_atom_queue;

    /**
     * m_container is the class which hold the actual cell lists
     *
     * it also manages iterating through the pairs of atoms
     */
    CellListsContainer<m_ndim> m_container;
public:
    ~CellLists() {}

    /**
     * constructor
     *
     * ncellx_scale scales the number of cells.  The number of cells in each
     * direction is computed from ncellx_scale * box_length / rcut
     */
    CellLists(
        std::shared_ptr<distance_policy> const & dist,
        pele::Array<double> const & boxv,
        const double rcut,
        const double ncellx_scale=1.0,
        const bool balance_omp=true);

    /**
     * return the class which loops over the atom pairs with a callback function
     */
    template <class callback_class>
    inline CellListsLoop<callback_class, m_ndim> get_atom_pair_looper(callback_class & callback) const
    {
        return CellListsLoop<callback_class, m_ndim>(callback, m_container);
    }

    /**
     * return the total number of cells
     */
    size_t get_nr_cells() const { return m_lattice_tool.m_ncells; }

    /**
     * return the number of cells in the x direction
     */
    size_t get_nr_cellsx() const { return m_lattice_tool.m_ncells_vec[0]; }

    pele::Array<size_t> get_order() const
    {
        auto order = pele::Array<size_t>(m_container.m_atoms.size());

        size_t ind = 0;
        for (size_t isubdom = 0; isubdom < m_lattice_tool.m_nsubdoms; ++isubdom) {
            for (size_t icell = 0; icell < m_lattice_tool.m_subdom_ncells[isubdom]; ++icell) {
                for (auto cell_iter = m_container.getIterator(icell, isubdom);
                     *cell_iter != CELL_END;
                     ++cell_iter) {
                    order[ind] = *cell_iter;
                    ++ind;
                }
            }
        }
        return order;
    }

    /**
     * update the cell list iterator with new coordinates
     */
    void update(pele::Array<double> const & coords);

protected:
    void print_warnings(const size_t natoms);
    void build_cell_neighbors_list();
    void reset_container(pele::Array<double> const & coords);
    void update_container(pele::Array<double> const & coords);
private:
    static Array<size_t> get_ncells_vec(Array<double> const & boxv, const double rcut,
                                        const double ncellx_scale, const bool balance_omp);
};

template<typename distance_policy>
CellLists<distance_policy>::CellLists(
        std::shared_ptr<distance_policy> const & dist,
        pele::Array<double> const & boxv,
        const double rcut,
        const double ncellx_scale,
        const bool balance_omp)
    : m_initialized(false),
      m_lattice_tool(dist, boxv, rcut, get_ncells_vec(boxv, rcut, ncellx_scale, balance_omp)),
      m_container(m_lattice_tool.m_subdom_ncells),
      add_atom_queue(m_lattice_tool.m_nsubdoms)
{
    build_cell_neighbors_list();

    if (m_ndim < 2) {
        throw std::runtime_error("CellLists::CellLists: Cell lists only work with more than 1 dimension (due to the split of subdomains in y-dimension)");
    }
    if (boxv.size() != m_ndim) {
        throw std::runtime_error("CellLists::CellLists: distance policy boxv and cell list boxv differ in size");
    }
    if (*std::min_element(boxv.begin(), boxv.end()) < rcut) {
        throw std::runtime_error("CellLists::CellLists: illegal rcut");
    }
    for (size_t i = 1; i < boxv.size(); ++i) {
        if (boxv[i] < 0) {
            throw std::runtime_error("CellLists::CellLists: illegal input: boxvector");
        }
    }
    if (ncellx_scale < 0) {
        throw std::runtime_error("CellLists::CellLists: illegal input: ncellx_scale");
    }
    #ifdef _OPENMP
    if (std::floor(ncellx_scale) != ncellx_scale && ncellx_scale > 1 && omp_get_max_threads() > 1) {
        throw std::runtime_error("CellLists::CellLists: Non-integer values > 1 of "
                                 "ncellx_scale can break the parallelization!");
    }
    #endif
    size_t ncell_min = *std::min_element(m_lattice_tool.m_ncells_vec.begin(), m_lattice_tool.m_ncells_vec.end());
    if (ncell_min < 5) {
        // If there are only a few cells in any direction then it doesn't make sense to use cell lists
        // because so many cells will be neighbors with each other.
        // It would be better to use simple loops over atom pairs.
        std::cout << "CellLists: efficiency warning: there are not many cells ("<<ncell_min<<") in at least one direction.\n";
    }
    if (m_lattice_tool.m_rcut > 0.5 * *std::min_element(m_lattice_tool.m_boxvec.begin(), m_lattice_tool.m_boxvec.end())) {
        // an atom can interact with more than just the nearest image of it's neighbor
        std::cerr << "CellLists: warning: rcut > half the box length.  This might cause errors with periodic boundaries.\n";
    }
}

template<typename distance_policy>
Array<size_t> CellLists<distance_policy>::get_ncells_vec(Array<double> const & boxv, const double rcut, const double ncellx_scale, const bool balance_omp)
{
    pele::Array<size_t> res(boxv.size());
    for (size_t i = 0; i < res.size(); ++i) {
        res[i] = std::max<size_t>(1, ncellx_scale * boxv[i] / rcut) ;
    }
    #ifdef _OPENMP
    if(omp_get_max_threads() > res[1]) {
        throw std::runtime_error("More threads than cells in y-direction. "
                                 "Reduce the number of threads "
                                 "(environment variable OMP_NUM_THREADS)!");
    }
    if(balance_omp) {
        res[1] = (res[1] / omp_get_max_threads()) * omp_get_max_threads();
    }
    #endif
    return res;
}

template<typename distance_policy>
void CellLists<distance_policy>::print_warnings(const size_t natoms)
{
    if (m_lattice_tool.m_ncells > natoms) {
        // It would be more efficient (I think) to reduce the number of cells.
        std::cout << "CellLists: efficiency warning: the number of cells ("<<m_lattice_tool.m_ncells<<")"<<
                " is greater than the number of atoms ("<<natoms<<").\n";
    }
}

/**
 * update or re-build cell lists
 */
template <typename distance_policy>
void CellLists<distance_policy>::update(pele::Array<double> const & coords)
{
    if (m_initialized) {
        update_container(coords);
    } else {
        m_initialized = true;
        size_t natoms = coords.size() / m_ndim;
        print_warnings(natoms);
        reset_container(coords);
    }
}

/**
 * build the list of neighboring cells.
 */
template <typename distance_policy>
void CellLists<distance_policy>::build_cell_neighbors_list()
{
    m_lattice_tool.find_neighbor_pairs(
        m_container.m_cell_neighbor_pairs_inner,
        m_container.m_cell_neighbor_pairs_boundary,
        m_container.m_cells);
}

/**
 * remove all atoms and re-add them to the right cells according to the new coordinates
 */
template <typename distance_policy>
void CellLists<distance_policy>::reset_container(pele::Array<double> const & coords)
{
    size_t natoms = coords.size() / m_ndim;
    m_container.m_atoms = std::vector<long>(natoms);
    m_container.reset();
    for(long iatom = 0; iatom < natoms; ++iatom) {
        const double * const x = coords.data() + m_ndim * iatom;
        size_t icell, isubdom;
        m_lattice_tool.position_to_local_ind(x, icell, isubdom);
        m_container.add_atom_to_cell(iatom, icell, isubdom);
    }
}

/**
 * re-calculate the cells for all atoms according to new coordinates
 */
template <typename distance_policy>
void CellLists<distance_policy>::update_container(pele::Array<double> const & coords)
{
    #ifdef _OPENMP
    #pragma omp parallel
    {
        size_t isubdom = omp_get_thread_num();
        for (size_t icell = 0; icell < m_lattice_tool.m_subdom_ncells[isubdom]; ++icell) {
            long parent_iatom = CELL_END;
            AtomInCellIterator<m_ndim> cell_iterator = m_container.getIterator(icell, isubdom);
            while (*cell_iterator != CELL_END) {
                long iatom = *cell_iterator;
                const double * const new_x = coords.data() + m_ndim * iatom;
                size_t new_cell, new_subdom;
                m_lattice_tool.position_to_local_ind(new_x, new_cell, new_subdom);
                if (new_cell != icell || new_subdom != isubdom) {
                    ++cell_iterator;
                    m_container.remove_atom_from_cell(parent_iatom, iatom, icell, isubdom);
                    if (isubdom == new_subdom) {
                        m_container.add_atom_to_cell(iatom, new_cell, isubdom);
                    } else {
                        std::array<long, 2> add_info = {iatom, new_cell};
                        add_atom_queue[new_subdom].push(add_info);
                    }
                } else {
                    ++cell_iterator;
                    parent_iatom = iatom;
                }
            }
        }

        #pragma omp barrier

        while (!add_atom_queue[isubdom].empty()) {
            auto add_info = add_atom_queue[isubdom].front();
            add_atom_queue[isubdom].pop();
            m_container.add_atom_to_cell(add_info[0], add_info[1], isubdom);
        }
    }
    #else
    for (size_t isubdom = 0; isubdom < m_lattice_tool.m_nsubdoms; ++isubdom) {
        for (size_t icell = 0; icell < m_lattice_tool.m_subdom_ncells[isubdom]; ++icell) {
            long parent_iatom = CELL_END;
            AtomInCellIterator<m_ndim> cell_iterator = m_container.getIterator(icell, isubdom);
            while (*cell_iterator != CELL_END) {
                long iatom = *cell_iterator;
                const double * const new_x = coords.data() + m_ndim * iatom;
                size_t new_cell, new_subdom;
                m_lattice_tool.position_to_local_ind(new_x, new_cell, new_subdom);
                if (new_cell != icell || new_subdom != isubdom) {
                    ++cell_iterator;
                    m_container.remove_atom_from_cell(parent_iatom, iatom, icell, isubdom);
                    m_container.add_atom_to_cell(iatom, new_cell, new_subdom);
                } else {
                    ++cell_iterator;
                    parent_iatom = iatom;
                }
            }
        }
    }
    #endif
}

} // namespace pele

#endif // #ifndef _PELE_CELL_LISTS_H_
