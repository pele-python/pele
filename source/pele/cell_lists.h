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

}

namespace pele{

template <size_t ndim>
struct AtomPosition {
    long atom_index;
    pele::VecN<ndim> x;

    AtomPosition()
        : atom_index(CELL_END)
    {}
    AtomPosition(long atom_index_, double const *x_)
        : atom_index(atom_index_),
          x(x_, x_+ndim)
    {}
    AtomPosition(AtomPosition<ndim> const & v)
        : atom_index(v.atom_index),
          x(v.x)
    {}
    inline bool operator!=(AtomPosition<ndim> const & v) const
    {
        return atom_index != v.atom_index;
    }
};
}

namespace {

/**
 * this iterator facilitates looping through the atoms in a cell.
 */
template<size_t ndim>
class AtomInCellIterator {
private:
    std::vector<pele::AtomPosition<ndim> > const * m_ll;
    pele::AtomPosition<ndim> m_current_atom;
public:
    AtomInCellIterator(std::vector<pele::AtomPosition<ndim> > const & ll,
            pele::AtomPosition<ndim> const & first_atom)
        : m_ll(&ll),
          m_current_atom(first_atom)
    {}

    AtomInCellIterator()
        : m_ll(NULL)
    {}

    /**
     * return the current atom
     */
    inline pele::AtomPosition<ndim> const & operator*() const
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
        m_current_atom = (*m_ll)[m_current_atom.atom_index];
    }

    inline bool operator!=(AtomInCellIterator<ndim> const & iter) const
    {
        return m_current_atom != iter.m_current_atom;
    }

};


/**
 * container for the cell lists
 */
template<size_t ndim>
class CellListsContainer {
public:
    /**
     * m_ll is an array of linked atoms.
     *
     * m_ll[atom_i.atom_index] is the index of the next atom in the same cell as atom_i.
     * if the index of m_ll[atom_i.atom_index] is CELL_END then there are no more atoms in this cell.
     */
    std::vector<pele::AtomPosition<ndim> > m_ll;

    /**
     * m_hoc is a head of chain list.
     *
     * m_hoc[icell] is the first atom in cell icell.  This is
     * used in conjuction with m_ll
     */
    std::vector<pele::AtomPosition<ndim> > m_hoc;

    /** vector of pairs of neighboring cells */
    std::vector<std::pair<size_t, size_t> > m_cell_neighbor_pairs;

    CellListsContainer(size_t ncells)
        : m_hoc(ncells)
    {}

    /**
     * clear all data from the lists
     */
    void clear()
    {
        for (auto & v : m_hoc) {
            v.atom_index = CELL_END;
        }
    }

    /**
     * add an atom to a cell
     */
    inline void add_atom_to_cell(size_t iatom, size_t icell, double const * atom_position)
    {
        m_ll[iatom] = m_hoc[icell];
        m_hoc[icell] = pele::AtomPosition<ndim>(iatom, atom_position);
//        std::cout << icell << " adding atom " << m_hoc[icell].atom_index << " " << m_hoc[icell].x << std::endl;
    }

    /**
     * set the size of the m_ll array
     *
     * this must be called before assigning the atoms to cells
     */
    inline void set_natoms(size_t natoms)
    {
        m_ll.resize(natoms);
    }

    typedef AtomInCellIterator<ndim> const_iterator;
    /**
     * return an iterator over the atoms in cell `icell`
     */
    AtomInCellIterator<ndim> begin(size_t icell) const
    {
        return AtomInCellIterator<ndim>(m_ll, m_hoc[icell]);
    }
    /**
     * return the end iterator for the atoms in a cell
     */
    AtomInCellIterator<ndim> end() const
    {
        return AtomInCellIterator<ndim>(m_ll, pele::AtomPosition<ndim>());
    }
};


} // end anonymous namespace


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

public:
    CellListsLoop(visitor_t & visitor, CellListsContainer<ndim> const & container)
        : m_visitor(visitor),
          m_container(container)
    {}

    void loop_through_atom_pairs()
    {
        typename CellListsContainer<ndim>::const_iterator iiter, jiter, iend, jend;
        iend = m_container.end();
        for (auto const & ijpair : m_container.m_cell_neighbor_pairs) {
            const size_t icell = ijpair.first;
            const size_t jcell = ijpair.second;
            // do double loop through atoms, avoiding duplicate pairs
            for (iiter = m_container.begin(icell); iiter != iend; ++iiter) {
                // if icell==jcell we need to avoid duplicate atom pairs
                jend = (icell == jcell) ? iiter : m_container.end();
                for (jiter = m_container.begin(jcell); jiter != jend; ++jiter) {
                    m_visitor.insert_atom_pair(*iiter, *jiter);
                }
            }
        }
    }
};


template<typename distance_policy>
class LatticeNeighbors {
public:
    static const size_t ndim = distance_policy::_ndim;
    std::shared_ptr<distance_policy> m_dist; // the distance function

    typedef VecN<ndim> cell_vec_t;
    pele::Array<double> m_boxvec;
    double m_rcut;
    cell_vec_t m_ncells_vec; // the number of cells in each dimension
    pele::VecN<ndim> m_rcell_vec; // the cell length in each dimension
    size_t m_ncells;

    LatticeNeighbors(std::shared_ptr<distance_policy> dist,
            pele::Array<double> boxvec,
            double rcut,
            pele::Array<size_t> ncells_vec)
        : m_dist(dist),
          m_boxvec(boxvec.copy()),
          m_rcut(rcut),
          m_ncells_vec(ncells_vec.begin(), ncells_vec.end())
    {
        m_ncells = m_ncells_vec.prod();
        for (size_t idim = 0; idim < ndim; ++idim) {
            m_rcell_vec[idim] = m_boxvec[idim] / m_ncells_vec[idim];
        }
    }

    /**
     * apply periodic boundary conditions to a cell vector
     */
    inline void cell_vec_apply_periodic(cell_vec_t & v) const
    {
        // note: the speed of this function is important because it is called
        // once per atom each time get_energy is called.
        for (size_t idim = 0; idim < ndim; ++idim) {
            v[idim] -= std::floor(v[idim] / m_ncells_vec[idim]) * m_ncells_vec[idim];
        }
    }

    /**
     * convert a cell vector to the cell index
     */
    size_t to_index(cell_vec_t v) const
    {
        // note: the speed of this function is important because it is called
        // once per atom each time get_energy is called.
        cell_vec_apply_periodic(v);
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
    cell_vec_t to_cell_vec(size_t icell) const
    {
        cell_vec_t v;
        size_t cum = m_ncells;
        for (long idim = ndim-1; idim >= 0 ; --idim) {
            cum /= m_ncells_vec[idim];
            v[idim] = size_t (icell / cum);
            icell -= v[idim] * cum;
        }
        return v;
    }

    /**
     * convert a cell vector to a positional vector in real space
     */
    VecN<ndim> to_position(cell_vec_t v) const
    {
        VecN<ndim> x;
        cell_vec_apply_periodic(v);
        for (size_t idim = 0; idim < ndim; ++idim) {
            x[idim] = m_rcell_vec[idim] * v[idim];
        }
        return x;
    }

    /**
     * return the index of the cell that contains position x
     */
    size_t position_to_cell_index(double const * const x) const
    {
        // note: the speed of this function is important because it is called
        // once per atom each time get_energy is called.
        cell_vec_t cell_vec;
        for(size_t idim = 0; idim < ndim; ++idim) {
            cell_vec[idim] = std::floor(m_ncells_vec[idim] * (x[idim]) / m_boxvec[idim]);
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
    size_t find_neighbors(size_t idim, cell_vec_t v0,
            std::vector<size_t> & neighbors,
            cell_vec_t const & vorigin
            ) const
    {
        if (idim == ndim) {
            double rmin = minimum_distance(v0, vorigin);
            if (rmin <= m_rcut) {
                neighbors.push_back(to_index(v0));
//                std::cout << "  adding neighbor    "<< v0 << " rmin " << rmin << "\n";
                return 1;
            }
//            std::cout << "  rejecting neighbor "<< v0 << " rmin " << rmin << "\n";
            return 0;
        }
//        std::cout << "find neighbors " << idim << " " << v0 << "\n";
        size_t nfound = 0;
        size_t n;
        auto v = v0;
        long max_negative = (m_ncells_vec[idim] - 1) / 2;
        long max_positive = m_ncells_vec[idim] / 2;
        // step in the negative direction
        long offset = 0;
        while (offset <= max_negative) {
            v[idim] = v0[idim] - offset;
            n = find_neighbors(idim+1, v, neighbors, vorigin);
            if (n == 0) break;
            nfound += n;
            ++offset;
        }
        // step in the positive direction
        offset = 1;
        while (offset <= max_positive) {
            v[idim] = v0[idim] + offset;
            n = find_neighbors(idim+1, v, neighbors, vorigin);
            if (n == 0) break;
            nfound += n;
            ++offset;
        }

        return nfound;
    }

    /**
     * return a vector of all the neighbors of icell (including icell itself)
     */
    std::vector<size_t> find_all_neighbors(size_t icell) const
    {
        auto vcell = to_cell_vec(icell);

        std::vector<size_t> neighbors;
        find_neighbors(0, vcell, neighbors, vcell);
//        std::cout << pele::Array<size_t>(neighbors) << std::endl;
        return neighbors;
    }

    /**
     * return a list of all pairs of neighboring cells
     *
     * The algorithm could be improved.  If the distance function is periodic
     * then each cell has the same structure of neighbors.  We could just keep
     * a list of cell vector offsets and apply these to each cell to find the
     * list of neighbor pairs.
     */
    void find_neighbor_pairs(std::vector<std::pair<size_t, size_t> > & cell_neighbors) const
    {
        cell_neighbors.reserve(m_ncells * std::pow(3, ndim));
        for (size_t icell = 0; icell < m_ncells; ++icell) {
            auto neighbors = find_all_neighbors(icell);
            for (size_t jcell : neighbors) {
                if (jcell >= icell) { // avoid duplicates
                    cell_neighbors.push_back(std::pair<size_t, size_t>(icell, jcell));
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

    pele::Array<double> m_coords; // the coordinates array
    size_t m_natoms; // the number of atoms
    bool m_initialised; // flag for whether the class has been initialized
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
     * direction is computed from ncellx_scale * box_lenth / rcut
     */
    CellLists(
        std::shared_ptr<distance_policy> dist,
        pele::Array<double> const boxv, const double rcut,
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
     * reset the cell list iterator with a new coordinates array
     */
    void reset(pele::Array<double> coords);

protected:
    void setup(Array<double> coords);
    void build_cell_neighbors_list();
    void build_linked_lists();
private:
    static Array<size_t> get_ncells_vec(const Array<double> boxv, const double rcut, const double ncellx_scale);
};



template<typename distance_policy>
CellLists<distance_policy>::CellLists(
        std::shared_ptr<distance_policy> dist,
        pele::Array<double> const boxv, const double rcut,
        const double ncellx_scale)
    : m_natoms(0),
      m_initialised(false),
      m_lattice_tool(dist, boxv, rcut, get_ncells_vec(boxv, rcut, ncellx_scale)),
      m_container(m_lattice_tool.m_ncells)
{
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

//    std::cout << "total number of cells " << m_ncells << std::endl;
}

template<typename distance_policy>
Array<size_t> CellLists<distance_policy>::get_ncells_vec(const Array<double> boxv, const double rcut, const double ncellx_scale)
{
    pele::Array<size_t> res(boxv.size());
    for (size_t i = 0; i < res.size(); ++i) {
        res[i] = std::max<size_t>(1, ncellx_scale * boxv[i] / rcut) ;
    }
    return res;
}

template<typename distance_policy>
void CellLists<distance_policy>::setup(Array<double> coords)
{
    m_coords = coords.copy();
    m_natoms = coords.size() / m_ndim;
    m_container.set_natoms(m_natoms);
    if (coords.size() != m_ndim * m_natoms) {
        throw std::runtime_error("CellLists::setup: illegal coords.size() not divisible by m_ndim");
    }

    build_cell_neighbors_list();
    m_initialised = true;

    // print messages if any of the parameters seem bad
    size_t ncell_min = *std::min_element(m_lattice_tool.m_ncells_vec.begin(), m_lattice_tool.m_ncells_vec.end());
    if (ncell_min < 5) {
        // If there are only a few cells in any direction then it doesn't make sense to use cell lists
        // because so many cells will be neighbors with each other.
        // It would be better to use simple loops over atom pairs.
        std::cout << "CellLists: efficiency warning: there are not many cells ("<<ncell_min<<") in at least one direction.\n";
    }
    if (m_lattice_tool.m_ncells > m_natoms) {
        // It would be more efficient (I think) to reduce the number of cells.
        std::cout << "CellLists: efficiency warning: the number of cells ("<<m_lattice_tool.m_ncells<<")"<<
                " is greater than the number of atoms ("<<m_natoms<<").\n";
    }
    if (m_lattice_tool.m_rcut > 0.5 * *std::min_element(m_lattice_tool.m_boxvec.begin(), m_lattice_tool.m_boxvec.end())) {
        // an atom can interact with more than just the nearest image of it's neighbor
        std::cerr << "CellLists: warning: rcut > half the box length.  This might cause errors with periodic boundaries.\n";
    }
}

/**
 * re-build linked lists
 * Algorithm 37 page 552 Understanding Molecular Simulation 2nd ed.
 * start by setting head of chain (hoc of size ncells) to -1 (meaning end of chain)
 * then update linked list so that atom i points to the next atom in the chain,
 * obviously this starts from -1 if it is the only element in the chain. If the next
 * atom i is in the same cell, then the hoc for that cell is set to be i
 * and the linked list at position i will point to the index of the previous atom.
 * This is done iteratively for all atoms.
 */
template <typename distance_policy>
void CellLists<distance_policy>::reset(pele::Array<double> coords)
{
    if (! m_initialised) {
        setup(coords);
    }

    m_coords.assign(coords);
//    if (periodic_policy_check<distance_policy>::is_periodic) {
//        // distance policy is periodic: put particles "back in box" first
//        periodic_distance<m_ndim>(m_lattice_tool.m_boxvec).put_in_box(m_coords);
//    }
//    else {
//        // distance policy is not periodic: check that particles are inside box
//        auto boxvec = m_lattice_tool.m_boxvec;
//        for (size_t i = 0; i < m_coords.size(); ++i) {
//            if (m_coords[i] < -0.5 * boxvec[0] || m_coords[i] > 0.5 * boxvec[0]) {
//                std::cout << "m_coords[i]: " << m_coords[i] << "\n";
//                std::cout << "0.5 * boxvec[0]: " << 0.5 * boxvec[0] << std::endl;
//                throw std::runtime_error("CellLists::reset: coords are incompatible with boxvector");
//            }
//        }
//    }
    build_linked_lists();
}

/**
 * build the list of neighboring cells.
 */
template <typename distance_policy>
void CellLists<distance_policy>::build_cell_neighbors_list()
{
    m_lattice_tool.find_neighbor_pairs(m_container.m_cell_neighbor_pairs);
}

/**
 * determine which cell each atom is in and populate the arrays hoc and ll
 */
template <typename distance_policy>
void CellLists<distance_policy>::build_linked_lists()
{
    m_container.clear();
    for(size_t iatom = 0; iatom < m_natoms; ++iatom) {
        double const * const x = m_coords.data() + m_ndim * iatom;
        size_t icell = m_lattice_tool.position_to_cell_index(x);
        m_container.add_atom_to_cell(iatom, icell, x);
    }
}

} // namespace pele





#endif // #ifndef _PELE_CELL_LISTS_H_
