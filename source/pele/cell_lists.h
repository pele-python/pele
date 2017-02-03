#ifndef _PELE_CELL_LISTS_H_
#define _PELE_CELL_LISTS_H_

#include <iostream>
#include <memory>
#include <exception>
#include <cassert>
#include <vector>

#include "base_potential.h"
#include "array.h"
#include "distance.h"
#include "vecn.h"

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
 * container for the cell lists
 */
template<size_t ndim>
class CellListsContainer {
public:
    /**
     * m_cell_atoms is a vector with vectors of atom indices.
     *
     * Each element (std::vector) contains the indices of all atoms in the corresponding cell.
     */
    std::vector< std::vector<size_t> > m_cell_atoms;

    /** vector of pairs of neighboring cells */
    std::vector<std::pair<const std::vector<size_t>*, const std::vector<size_t>*> > m_cell_neighbor_pairs;

    CellListsContainer(const size_t ncells)
        : m_cell_atoms(ncells)
    {
    }

    /**
     * clear all data from the vectors
     */
    void clear()
    {
        for (auto & v : m_cell_atoms) {
            v.clear();
        }
    }

    /**
     * add an atom to a cell
     */
    inline void add_atom_to_cell(const size_t iatom, const size_t icell)
    {
        m_cell_atoms[icell].push_back(iatom);
    }

    /**
     * remove an atom from a cell
     */
    inline void remove_atom_from_cell(const size_t iatom, const size_t icell)
    {
        auto index = std::find(m_cell_atoms[icell].begin(), m_cell_atoms[icell].end(), iatom);

        if (index != m_cell_atoms[icell].end()) {
            std::swap(*index, m_cell_atoms[icell].back());
        }
        m_cell_atoms[icell].pop_back();
    }

    /**
     * remove an atom from a cell using its index within that cell
     */
    inline void remove_atom_from_cell_by_index(const size_t index, const size_t icell)
    {
        std::swap(m_cell_atoms[icell][index], m_cell_atoms[icell].back());
        m_cell_atoms[icell].pop_back();
    }
};

}

namespace pele{

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

public:
    CellListsLoop(visitor_t & visitor, CellListsContainer<ndim> const & container)
        : m_visitor(visitor),
          m_container(container)
    {}

    void loop_through_atom_pairs(Array<double> const & coords)
    {
        for (auto const & ijpair : m_container.m_cell_neighbor_pairs) {
            const std::vector<size_t>* icell = ijpair.first;
            const std::vector<size_t>* jcell = ijpair.second;
            // do double loop through atoms, avoiding duplicate pairs
            for (auto iatom = icell->begin(); iatom != icell->end(); ++iatom) {
                // if icell==jcell we need to avoid duplicate atom pairs
                auto jend = (icell == jcell) ? iatom : jcell->end();
                for (auto jatom = jcell->begin(); jatom != jend; ++jatom) {
                    m_visitor.insert_atom_pair(coords, *iatom, *jatom);
                }
            }
        }
    }
};

template<typename distance_policy>
class LatticeNeighbors {
public:
    static const size_t ndim = distance_policy::_ndim;
    const std::shared_ptr<distance_policy> m_dist; // the distance function

    typedef VecN<ndim, size_t> cell_vec_t;
    pele::VecN<ndim> m_boxvec;
    pele::VecN<ndim> m_inv_boxvec; // inverse of boxvec
    double m_rcut;
    cell_vec_t m_ncells_vec; // the number of cells in each dimension
    pele::VecN<ndim> m_rcell_vec; // the cell length in each dimension
    size_t m_ncells;

    LatticeNeighbors(std::shared_ptr<distance_policy> const & dist,
            pele::Array<double> const & boxvec,
            const double rcut,
            pele::Array<size_t> const & ncells_vec)
        : m_dist(dist),
          m_boxvec(boxvec),
          m_rcut(rcut),
          m_ncells_vec(ncells_vec.begin(), ncells_vec.end()),
          m_inv_boxvec(ndim)
    {
        m_ncells = m_ncells_vec.prod();
        for (size_t idim = 0; idim < ndim; ++idim) {
            m_inv_boxvec[idim] = 1 / m_boxvec[idim];
            m_rcell_vec[idim] = m_boxvec[idim] / m_ncells_vec[idim];
        }
    }

    /**
     * convert a cell vector to the cell index
     * The cell vector needs to be within range [0, m_ncells_vec]
     */
    size_t to_index(cell_vec_t const & v) const
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
     * convert a cell index to a cell vector
     */
    cell_vec_t to_cell_vec(const size_t icell) const
    {
        cell_vec_t v;
        size_t remaining_icell = icell;
        size_t fraction;
        for (long idim = 0; idim < ndim - 1; idim++) {
            fraction = (size_t) remaining_icell / m_ncells_vec[idim];
            // should be optimized by compiler to use remainder of previous division
            v[idim] = remaining_icell % m_ncells_vec[idim];
            remaining_icell = fraction;
        }
        v[ndim - 1] = remaining_icell;
        return v;
    }

    /**
     * convert a cell vector to a positional vector in real space
     */
    VecN<ndim> to_position(cell_vec_t const & v) const
    {
        VecN<ndim> x;
        for (size_t idim = 0; idim < ndim; ++idim) {
            x[idim] = m_rcell_vec[idim] * v[idim];
        }
        return x;
    }

    /**
     * return the index of the cell that contains position x
     */
    size_t position_to_cell_index(double * const x) const
    {
        // note: the speed of this function is important because it is called
        // once per atom each time get_energy is called.
        cell_vec_t cell_vec;
        for(size_t idim = 0; idim < ndim; ++idim) {
            if (x[idim] <= -0.5 * m_boxvec[idim] || x[idim] >= 0.5 * m_boxvec[idim]) {
                m_dist->put_atom_in_box(x);
            }
            cell_vec[idim] = std::floor(m_ncells_vec[idim] * (x[idim] * m_inv_boxvec[idim] + 0.5));
        }
        return to_index(cell_vec);
    }

    /**
     * return the minimum corner to corner distance between two cells
     */
    double minimum_distance(cell_vec_t const & v1, cell_vec_t const & v2) const
    {
        // copy them so we don't accidentally change them
        auto lower_left1 = to_position(v1);
        auto lower_left2 = to_position(v2);
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
     */
    void find_neighbors(const size_t idim, cell_vec_t const & v0,
            std::vector<size_t> & neighbors,
            cell_vec_t const & vorigin
            ) const
    {
        if (idim == ndim) {
            double rmin = minimum_distance(v0, vorigin);
            if (rmin <= m_rcut) {
                neighbors.push_back(to_index(v0));
            }
        } else {
            auto v = v0;
            for(v[idim] = 0; v[idim] < m_ncells_vec[idim]; v[idim]++) {
                find_neighbors(idim+1, v, neighbors, vorigin);
            }
        }
    }

    /**
     * return a vector of all the neighbors of icell (including icell itself)
     */
    std::vector<size_t> find_all_neighbors(const size_t icell) const
    {
        auto vcell = to_cell_vec(icell);

        std::vector<size_t> neighbors;
        find_neighbors(0, vcell, neighbors, vcell);
        return neighbors;
    }

    /**
     * return a list of all pairs of neighboring cells
     *
     * The algorithm could be improved.  If the distance function is periodic
     * then each cell has the same structure of neighbors.  We could just keep
     * a list of cell vector offsets and apply these to each cell to find the
     * list of neighbor pairs.
     * This would not be true for Lees-Edwards boundary conditions, though.
     * Since this algorithm is not performance-critical, I would argue against
     * such a specialization.
     */
    void find_neighbor_pairs(
        std::vector< std::pair< const std::vector<size_t>*, const std::vector<size_t>* > > & cell_neighbors,
        std::vector< std::vector<size_t> > const & cells) const
    {
        cell_neighbors.reserve(m_ncells * std::pow(3, ndim));
        for (size_t icell = 0; icell < m_ncells; ++icell) {
            auto neighbors = find_all_neighbors(icell);
            for (size_t jcell : neighbors) {
                if (jcell >= icell) { // avoid duplicates
                    cell_neighbors.push_back(
                        std::pair<const std::vector<size_t>*, const std::vector<size_t>*>(
                        &cells[icell], &cells[jcell]));
                }
            }
        }
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
    typedef CellListsContainer<m_ndim> container_type;
protected:

    size_t m_natoms; // the number of atoms
    bool m_initialized; // flag for whether the cell lists have been initialized with coordinates
    pele::LatticeNeighbors<distance_policy> m_lattice_tool;

    /**
     * m_container is the class which hold the actual cell lists
     *
     * it also manages iterating through the pairs of atoms
     */
    container_type m_container;
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
        pele::Array<double> const & boxv, const double rcut,
        const double ncellx_scale=1.0);

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

    /**
     * update the cell list iterator with new coordinates
     */
    void update(pele::Array<double> & coords);

protected:
    void print_warnings(const size_t natoms);
    void build_cell_neighbors_list();
    void reset_container(pele::Array<double> & coords);
    void update_container(pele::Array<double> & coords);
private:
    static Array<size_t> get_ncells_vec(Array<double> const & boxv, const double rcut, const double ncellx_scale);
};

template<typename distance_policy>
CellLists<distance_policy>::CellLists(
        std::shared_ptr<distance_policy> const & dist,
        pele::Array<double> const & boxv, const double rcut,
        const double ncellx_scale)
    : m_initialized(false),
      m_lattice_tool(dist, boxv, rcut, get_ncells_vec(boxv, rcut, ncellx_scale)),
      m_container(m_lattice_tool.m_ncells)
{
    build_cell_neighbors_list();

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
        throw std::runtime_error("CellLists::CellLists: illegal input");
    }
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
Array<size_t> CellLists<distance_policy>::get_ncells_vec(Array<double> const & boxv, const double rcut, const double ncellx_scale)
{
    pele::Array<size_t> res(boxv.size());
    for (size_t i = 0; i < res.size(); ++i) {
        res[i] = std::max<size_t>(1, ncellx_scale * boxv[i] / rcut) ;
    }
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
 * re-build cell lists
 */
template <typename distance_policy>
void CellLists<distance_policy>::update(pele::Array<double> & coords)
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
        m_container.m_cell_neighbor_pairs,
        m_container.m_cell_atoms);
}

/**
 * remove all atoms and re-add them to the right cells according to the new coordinates
 */
template <typename distance_policy>
void CellLists<distance_policy>::reset_container(pele::Array<double> & coords)
{
    m_container.clear();
    size_t natoms = coords.size() / m_ndim;
    for(size_t iatom = 0; iatom < natoms; ++iatom) {
        double * const x = coords.data() + m_ndim * iatom;
        size_t icell = m_lattice_tool.position_to_cell_index(x);
        m_container.add_atom_to_cell(iatom, icell);
    }
}

/**
 * re-calculate the cells for all atoms according to new coordinates
 */
template <typename distance_policy>
void CellLists<distance_policy>::update_container(pele::Array<double> & coords)
{
    for (size_t icell = 0; icell < m_lattice_tool.m_ncells; ++icell) {
        size_t atom_nr = 0;
        while (atom_nr < m_container.m_cell_atoms[icell].size()) {
            size_t iatom = m_container.m_cell_atoms[icell][atom_nr];
            double * const new_x = coords.data() + m_ndim * iatom;
            size_t new_cell = m_lattice_tool.position_to_cell_index(new_x);
            if (new_cell != icell) {
                m_container.remove_atom_from_cell_by_index(atom_nr, icell);
                m_container.add_atom_to_cell(iatom, new_cell);
            } else {
                atom_nr++;
            }
        }
    }
}

} // namespace pele

#endif // #ifndef _PELE_CELL_LISTS_H_
