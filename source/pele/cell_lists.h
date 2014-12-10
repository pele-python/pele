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
 * this iterator facilitates looping through the atoms in a cell.
 *
 * It acts like an iterator in some ways, but it is not a real iterator.
 */
class AtomInCellIterator {
private:
    long const * const m_ll;
    long m_current_atom;
    long const m_end;
public:
    AtomInCellIterator(long const * const ll, size_t first_atom, long end=CELL_END)
        : m_ll(ll),
          m_current_atom(first_atom),
          m_end(end)
    {}
    AtomInCellIterator()
        : m_ll(NULL),
          m_current_atom(CELL_END),
          m_end(CELL_END)
    {}

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
        m_current_atom = m_ll[m_current_atom];
    }

    inline bool done() const
    {
        // If this were a real iterator, you would test if it's
        // done by comparing to some end() iterator.  I havn't quite
        // figured out how to do that in fast way.
        return m_current_atom == m_end;
    }
};


/**
 * container for the cell lists
 */
class CellListsContainer {
public:
    /**
     * m_ll is an array of linked atom indices.
     *
     * m_ll[atom_i] is the index of the next atom in the same cell as atom_i.
     * if m_ll[atom_i] is CELL_END then there are no more atoms in this cell.
     */
    pele::Array<long> m_ll;

    /**
     * m_hoc is a head of chain list.
     *
     * m_hoc[icell] is the index of the first atom in cell icell.  This is
     * used in conjuction with m_ll
     */
    pele::Array<long> m_hoc;

    /** vector of pairs of neighboring cells */
    std::vector<std::pair<size_t, size_t> > m_cell_neighbor_pairs;

    CellListsContainer(size_t ncells)
        : m_hoc(ncells, CELL_END)
    {}

    /**
     * clear all data from the lists
     */
    void clear()
    {
        m_hoc.assign(CELL_END);
    }

    /**
     * add an atom to a cell
     */
    inline void add_atom_to_cell(size_t iatom, size_t icell)
    {
        m_ll[iatom] = m_hoc[icell];
        m_hoc[icell] = iatom;
    }

    /**
     * set the size of the m_ll array
     */
    inline void set_natoms(size_t natoms)
    {
        m_ll = pele::Array<long>(natoms);
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
template <class visitor_t>
class CellListsLoop {
protected:
    visitor_t & m_visitor;
    pele::Array<long> const m_ll;
    pele::Array<long> const m_hoc;
    std::vector<std::pair<size_t, size_t> > const & m_cell_neighbor_pairs;

public:
    virtual ~CellListsLoop() {}
    
    CellListsLoop(visitor_t & visitor, CellListsContainer const & container)
        : m_visitor(visitor),
          m_ll(container.m_ll),
          m_hoc(container.m_hoc),
          m_cell_neighbor_pairs(container.m_cell_neighbor_pairs)
    {}


    void loop_through_atom_pairs()
    {
        for (auto const & ijpair : m_cell_neighbor_pairs) {
            const size_t icell = ijpair.first;
            const size_t jcell = ijpair.second;
            // do double loop through atoms, avoiding duplicate pairs
            for (auto iiter = AtomInCellIterator(m_ll.data(), m_hoc[icell]); !iiter.done(); ++iiter) {
                size_t const atomi = *iiter;
                // if icell==jcell we need to avoid duplicate atom pairs
                long const loop_end = (icell == jcell) ? atomi : CELL_END;
                for (auto jiter = AtomInCellIterator(m_ll.data(), m_hoc[jcell], loop_end); !jiter.done(); ++jiter) {
                    size_t const atomj = *jiter;
                    m_visitor.insert_atom_pair(atomi, atomj);
                }
            }
        }

    }

};

/**
 * Looping over atom pair similar to CellListsLoop.
 * The loop in loop_through_atom_pairs can terminate based on the result
 * of m_visitor.insert_atom_pair(atomi, atomj).
 */
template <class visitor_t>
class CellListsLoopBreak : public CellListsLoop<visitor_t> {
public:
    virtual ~CellListsLoopBreak() {}
    CellListsLoopBreak(visitor_t& visitor, CellListsContainer const& container)
        : CellListsLoop<visitor_t>(visitor, container)
    {}
    void loop_through_atom_pairs()
    {
        for (auto const & ijpair : CellListsLoop<visitor_t>::m_cell_neighbor_pairs) {
            const size_t icell = ijpair.first;
            const size_t jcell = ijpair.second;
            // do double loop through atoms, avoiding duplicate pairs
            for (auto iiter = AtomInCellIterator(CellListsLoop<visitor_t>::m_ll.data(), CellListsLoop<visitor_t>::m_hoc[icell]); !iiter.done(); ++iiter) {
                size_t const atomi = *iiter;
                // if icell==jcell we need to avoid duplicate atom pairs
                long const loop_end = (icell == jcell) ? atomi : CELL_END;
                for (auto jiter = AtomInCellIterator(CellListsLoop<visitor_t>::m_ll.data(), CellListsLoop<visitor_t>::m_hoc[jcell], loop_end); !jiter.done(); ++jiter) {
                    size_t const atomj = *jiter;
                    const bool break_loop = CellListsLoop<visitor_t>::m_visitor.insert_atom_pair(atomi, atomj);
                    if (break_loop) {
                        return;
                    }
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
    typedef CellListsContainer container_type;
protected:
    static const size_t m_ndim = distance_policy::_ndim;

    std::shared_ptr<distance_policy> m_dist; // the distance function
    pele::Array<double> m_coords; // the coordinates array
    size_t m_natoms; // the number of atoms
    const double m_rcut; // the potential cutoff
    bool m_initialised; // flag for whether the class has been initialized
    const pele::Array<double> m_boxv; // the array of box lengths
    const size_t m_ncellx; // the number of cells in the x direction
    const size_t m_ncells; // the total number of cells
    const double m_rcell; // the size of a cell

    /**
     * m_container is the class which hold the actual cell lists
     *
     * it also manages iterating through the pairs of atoms
     */
    container_type m_container;
    const double m_xmin;
    const double m_xmax;
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
    inline CellListsLoop<callback_class> get_atom_pair_looper(callback_class & callback) const
    {
        return CellListsLoop<callback_class>(callback, m_container);
    }

    /**
     * return the total number of cells
     */
    size_t get_nr_cells() const { return m_ncells; }

    /**
     * return the number of cells in the x direction
     */
    size_t get_nr_cellsx() const { return m_ncellx; }

    /**
     * return the number of unique atom pairs.
     *
     * These three functions are primarily used for debugging and testing
     */
    size_t get_nr_unique_pairs() const;
    size_t get_direct_nr_unique_pairs(const double max_distance, pele::Array<double> x) const;
    size_t get_maximum_nr_unique_pairs(pele::Array<double> x) const;

    /**
     * reset the cell list iterator with a new coordinates array
     */
    void reset(pele::Array<double> coords);

protected:
    void setup(Array<double> coords);
    void sanity_check();
    size_t atom2xbegin(const size_t atom_index) const { return m_ndim * atom_index; }
    template <class T> T loop_pow(const T x, int ex) const;
    size_t atom2cell(const size_t i);
    pele::Array<double> cell2coords(const size_t icell) const;
    bool cells_are_neighbors(const size_t icell, const size_t jcell) const;
    double get_minimum_corner_distance2(pele::Array<double> ic, pele::Array<double> jc) const;
    void build_cell_neighbors_list();
    void build_atom_neighbors_list_old();
    void build_atom_neighbors_list();
    void build_linked_lists();
};



template<typename distance_policy>
CellLists<distance_policy>::CellLists(
        std::shared_ptr<distance_policy> dist,
        pele::Array<double> const boxv, const double rcut,
        const double ncellx_scale)
    : m_dist(dist),
      m_natoms(0),
      m_rcut(rcut),
      m_initialised(false),
      m_boxv(boxv.copy()),
      m_ncellx(std::max<size_t>(1, (size_t)(ncellx_scale * m_boxv[0] / rcut))),     //no of cells in one dimension
      m_ncells(std::pow(m_ncellx, m_ndim)),                                                     //total no of cells
      m_rcell(m_boxv[0] / static_cast<double>(m_ncellx)),                                      //size of cell
      m_container(m_ncells),                                                                         //head of chain
      m_xmin(-0.5 * m_boxv[0]),
      m_xmax(0.5 * m_boxv[0])                                     
{
    if (m_boxv.size() != m_ndim) {
        throw std::runtime_error("CellLists::CellLists: distance policy boxv and cell list boxv differ in size");
    }
    if (*std::min_element(m_boxv.data(), m_boxv.data() + m_ndim) < rcut) {
        throw std::runtime_error("CellLists::CellLists: illegal rcut");
    }
    const double boxv_epsilon = 1e-10;
    for (size_t i = 1; i < boxv.size(); ++i) {
        if (fabs(boxv[0] - boxv[i]) > boxv_epsilon) {
            throw std::runtime_error("CellLists::CellLists: illegal input boxv is not for square box");
        }
        if (boxv[i] < 0) {
            throw std::runtime_error("CellLists::CellLists: illegal input: boxvector");
        }
    }
    if (m_ncellx == 0) {
        throw std::runtime_error("CellLists::CellLists: illegal lattice spacing");
    }
    if (ncellx_scale < 0) {
        throw std::runtime_error("CellLists::CellLists: illegal input");
    }

//    std::cout << "total number of cells " << m_ncells << std::endl;
}

// this is used only for tests.  It should be moved to the tests folder
class stupid_counter {
public:
    size_t count;
    stupid_counter() : count(0) {}
    void insert_atom_pair(size_t const atom_i, size_t const atom_j)
    { ++count; }
};

template<typename distance_policy>
size_t CellLists<distance_policy>::get_nr_unique_pairs() const
{
    // this function should be removed
    stupid_counter counter;
    auto looper = this->get_atom_pair_looper(counter);
    looper.loop_through_atom_pairs();
    return counter.count;
}


template<typename distance_policy>
size_t CellLists<distance_policy>::get_direct_nr_unique_pairs(const double max_distance, pele::Array<double> x) const
{
    size_t nr_unique_pairs = 0;
    const size_t natoms = x.size() / m_ndim;
    for (size_t i = 0; i < natoms; ++i) {
        for (size_t j = i + 1; j < natoms; ++j) {
            double rij[m_ndim];
            const double* xi = x.data() + atom2xbegin(i);
            const double* xj = x.data() + atom2xbegin(j);
            m_dist->get_rij(rij, xi, xj);
            double r2 = 0;
            for (size_t k = 0; k < m_ndim; ++k) {
                r2 += rij[k] * rij[k];
            }
            nr_unique_pairs += (r2 <= (max_distance * max_distance));
        }
    }
    return nr_unique_pairs;
}

template <typename distance_policy>
size_t CellLists<distance_policy>::get_maximum_nr_unique_pairs(pele::Array<double> x) const
{
    const size_t natoms = x.size() / m_ndim;
    return (natoms * (natoms - 1)) / 2;
}

template <typename distance_policy>
void CellLists<distance_policy>::setup(Array<double> coords)
{
    m_coords = coords.copy();
    m_natoms = coords.size() / m_ndim;
    m_container.set_natoms(m_natoms);
    if (coords.size() != m_ndim * m_natoms) {
        throw std::runtime_error("CellLists::setup: illegal coords.size() not divisible by m_ndim");
    }

//    m_atom_neighbor_list.reserve(m_natoms * (m_natoms - 1) / 2);
    build_cell_neighbors_list();
    m_initialised = true;
//    sanity_check();

    // print messages if any of the parameters seem bad
    if (m_ncellx < 5) {
        // If there are only a few cells in any direction then it doesn't make sense to use cell lists
        // because so many cells will be neighbors with each other.
        // It would be better to use simple loops over atom pairs.
        std::cout << "CellLists: efficiency warning: there are not many cells ("<<m_ncellx<<") in each direction.\n";
    }
    if (m_ncells > m_natoms) {
        // It would be more efficient (I think) to reduce the number of cells.
        std::cout << "CellLists: efficiency warning: the number of cells ("<<m_ncells<<")"<<
                " is greater than the number of atoms ("<<m_natoms<<").\n";
    }
    if (m_rcut > 0.5 * m_boxv[0]) {
        // an atom can interact with more than just the nearest image of it's neighbor
        std::cerr << "CellLists: warning: rcut > half the box length.  This might cause errors with periodic boundaries.\n";
    }
}

template <typename distance_policy>
void CellLists<distance_policy>::sanity_check()
{
    const size_t nr_unique_pairs_lists = get_nr_unique_pairs();
    const size_t nr_unique_pairs_direct = get_direct_nr_unique_pairs(m_rcut, m_coords);
    const size_t maximum_nr_unique_pairs = get_maximum_nr_unique_pairs(m_coords);
    //std::cout << "nr_unique_pairs_lists: " << nr_unique_pairs_lists << "\n";
    //std::cout << "nr_unique_pairs_direct: " << nr_unique_pairs_direct << "\n";
    //std::cout << "maximum_nr_unique_pairs: " << maximum_nr_unique_pairs << "\n";
    if (nr_unique_pairs_lists < nr_unique_pairs_direct) {
        std::cout << "nr_unique_pairs_lists: " << nr_unique_pairs_lists << "\n";
        std::cout << "nr_unique_pairs_direct: " << nr_unique_pairs_direct << "\n";
        std::cout << "maximum_nr_unique_pairs: " << maximum_nr_unique_pairs << "\n";
        throw std::runtime_error("CellLists::setup: sanity check failed: too few pairs");
    }
    if (nr_unique_pairs_lists > maximum_nr_unique_pairs) {
        throw std::runtime_error("CellLists::setup: sanity check failed: too many pairs");
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
    if (periodic_policy_check<distance_policy>::is_periodic) {
        // distance policy is periodic: put particles "back in box" first
        periodic_distance<m_ndim>(m_boxv).put_in_box(m_coords);
    }
    else {
        // distance policy is not periodic: check that particles are inside box
        for (size_t i = 0; i < m_coords.size(); ++i) {
            if (m_coords[i] < -0.5 * m_boxv[0] || m_coords[i] > 0.5 * m_boxv[0]) {
                std::cout << "m_coords[i]: " << m_coords[i] << "\n";
                std::cout << "0.5 * m_boxv[0]: " << 0.5 * m_boxv[0] << std::endl;
                throw std::runtime_error("CellLists::reset: coords are incompatible with boxvector");
            }
        }
    }
    build_linked_lists();
}

/**
 * return x to the power ex
 *
 * This is equivalent, but hopefully faster than std::pow(x, ex)
 */
template <typename distance_policy>
template <class T>
T CellLists<distance_policy>::loop_pow(const T x, int ex) const
{
    T result(1);
    while (--ex > -1) {
        result *= x;
    }
    return result;
}

/**
 * return the index of the cell that atom i is in
 *
 * this function assumes that particles have been already put in box
 */
template <typename distance_policy>
size_t CellLists<distance_policy>::atom2cell(const size_t i)
{
    assert(i < m_natoms);
    size_t icell = 0;
    for(size_t j = 0; j < m_ndim; ++j) {
        // j1 is the index for the coords array
        const size_t j1 = atom2xbegin(i) + j;
        assert(j1 < m_coords.size());
        double x = m_coords[j1];
        // min is needed in case x == m_rcell * m_ncellx
        const size_t icell_jpart = std::min<size_t>(m_ncellx - 1, static_cast<size_t>(((x - m_xmin) / (m_xmax - m_xmin)) * m_ncellx));
        assert(icell_jpart == icell_jpart);
        if (icell_jpart >= m_ncellx) {
            std::cout << "x: " << x << std::endl;
            std::cout << "m_rcell: " << m_rcell << std::endl;
            std::cout << "m_ndim: " << m_ndim << std::endl;
            std::cout << "m_ncellx: " << m_ncellx << std::endl;
            std::cout << "icell_jpart: " << icell_jpart << std::endl;
        }
        assert(icell_jpart < m_ncellx);
        //icell += icell_jpart * std::pow(m_ncellx, j);
        icell += icell_jpart * loop_pow(m_ncellx, j);
    }
    assert(icell < m_ncells);
    return icell;
}

/**
 * returns the coordinates to the corner of the lower left corner of cell icell
 *
 * lower-left means that the cartesian coordinates are smaller than all other corners.
 */
template <typename distance_policy>
pele::Array<double> CellLists<distance_policy>::cell2coords(const size_t icell) const
{
    pele::Array<double> cellcorner(m_ndim); //coordinate of cell bottom left corner
    std::vector<double> indexes(m_ndim, 0); //this array will store indexes, set to 0
    double index = 0;

    //don't change these loops to size_t or the conditions will not hold
    for(int i = m_ndim - 1; i >= 0; --i) {
        index = icell;
        for (int j = m_ndim - 1; j >= i; --j) {
            //index -= indexes[j] * std::pow(m_ncellx, j);
            index -= indexes[j] * loop_pow(m_ncellx, j);
        }
        //indexes[i] = floor(index / std::pow(m_ncellx, i));
        indexes[i] = floor(index / loop_pow(m_ncellx, i));
        cellcorner[i] = m_rcell * indexes[i];
    }

    return cellcorner.copy();
}

/**
 * return true if the cells are neighbors.
 *
 * The cells are considered neighbors if atoms in the cells could possibly
 * be closer than the cutoff distance
 */
template <typename distance_policy>
bool CellLists<distance_policy>::cells_are_neighbors(const size_t icell, const size_t jcell) const
{
    if (icell == jcell) {
        return true;
    }
    // Get "lower-left" corners.
    pele::Array<double> icell_coords = cell2coords(icell);
    pele::Array<double> jcell_coords = cell2coords(jcell);
    return get_minimum_corner_distance2(icell_coords, jcell_coords) <= m_rcut * m_rcut;
}

template <typename distance_policy>
double CellLists<distance_policy>::get_minimum_corner_distance2(pele::Array<double> lower_left1, pele::Array<double> lower_left2) const
{
    // copy them so we don't accidentally change them
    lower_left1 = lower_left1.copy();
    lower_left2 = lower_left2.copy();
    pele::VecN<m_ndim> ll1, ll2, dr;
    pele::VecN<m_ndim> minimum_distance; // the minimum possible distance in each direction
    for (size_t i = 0; i < m_ndim; ++i) {
        double min_dist = std::numeric_limits<double>::max();
        double dri;
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
        ll1[i] += m_rcell;
        m_dist->get_rij(dr.data(), ll1.data(), ll2.data());
        dri = std::abs(dr[i]);
        if (dri < min_dist) {
            min_dist = dri;
        }

        ll1 = lower_left1;
        ll2 = lower_left2;
        ll2[i] += m_rcell;
        m_dist->get_rij(dr.data(), ll1.data(), ll2.data());
        dri = std::abs(dr[i]);
        if (dri < min_dist) {
            min_dist = dri;
        }

        ll1 = lower_left1;
        ll2 = lower_left2;
        ll1[i] += m_rcell;
        ll2[i] += m_rcell;
        m_dist->get_rij(dr.data(), ll1.data(), ll2.data());
        dri = std::abs(dr[i]);
        if (dri < min_dist) {
            min_dist = dri;
        }

        minimum_distance[i] = min_dist;
    }
    double r2_min = dot<m_ndim> (minimum_distance, minimum_distance);
    return r2_min;
}

/**
 * build the list of neighboring cells.
 */
template <typename distance_policy>
void CellLists<distance_policy>::build_cell_neighbors_list()
{
    size_t max_n_neibs = 0;
    m_container.m_cell_neighbor_pairs.reserve(2 * m_ncells); // A lower end guess for the size
    for(size_t i = 0; i < m_ncells; ++i) {
        size_t nneibs = -0;
        for(size_t j = 0; j <= i; ++j) {
            if (cells_are_neighbors(i, j)) { //includes itself as a neighbor
                m_container.m_cell_neighbor_pairs.push_back(std::pair<size_t, size_t>(i, j));
                ++nneibs;
            }
        }
        if (nneibs > max_n_neibs) max_n_neibs = nneibs;
    }
    if (max_n_neibs > 0.5 * m_ncells) {
        // If each cell has many neighbors it would be better to just use a simple loop over atom pairs.
        // Alternatively you might think abour reducing rcut.
        std::cout << "CellLists: efficiency warning: the cells have very many neighbors ("
                <<max_n_neibs << ", with "<<m_ncells<<" cells total).\n";
    }

}

/**
 * determine which cell each atom is in and populate the arrays hoc and ll
 */
template <typename distance_policy>
void CellLists<distance_policy>::build_linked_lists()
{
    m_container.clear();
    for(size_t i = 0; i < m_natoms; ++i) {
        size_t icell = atom2cell(i);
        m_container.add_atom_to_cell(i, icell);
    }
}

} // namespace pele





#endif // #ifndef _PELE_CELL_LISTS_H_
