#ifndef PELE_NEIGHBOR_ITERATOR_H
#define PELE_NEIGHBOR_ITERATOR_H

#include <iostream>
#include <memory>
#include <exception>
#include <cassert>
#include <vector>

#include "base_potential.h"
#include "array.h"
#include "distance.h"

namespace pele {

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

/*cell list currently only work with box of equal side lengths
 * cell lists are currently not implemented for non cubic boxes:
 * in that case m_ncellx should be an array rather than a scalar and the definition
 * of ncells and rcell would change to array. This implies that in all the function these
 * would need to be replace with the correct array element. This adds room for errors
 * so in this first implementation we do not account for that scenario
 * */

template<typename distance_policy = periodic_distance<3> >
class CellIter{
public:
    typedef std::vector<std::pair<size_t, size_t> > container_type;
    typedef typename container_type::const_iterator const_iterator;
protected:
    static const size_t m_ndim = distance_policy::_ndim;
    std::shared_ptr<distance_policy> m_dist;
    pele::Array<double> m_coords;
    const size_t m_natoms;
    const double m_rcut;
    bool m_initialised;
    const pele::Array<double> m_boxv;
    const size_t m_ncellx;
    const size_t m_ncells;
    const double m_rcell;
    pele::Array<long int> m_hoc;
    pele::Array<long int> m_ll;
    std::vector<std::vector<size_t> > m_cell_neighbors;
    std::vector<std::pair<size_t, size_t> > m_atom_neighbor_list;
    const_iterator m_container_iterator;
    const double m_xmin;
    const double m_xmax;
public:
    ~CellIter() {}
    CellIter(pele::Array<double> const coords,
        std::shared_ptr<distance_policy> dist,
        pele::Array<double> const boxv, const double rcut,
        const double ncellx_scale=1.0);
    const_iterator begin() const { return m_atom_neighbor_list.begin(); }
    const_iterator end() const { return m_atom_neighbor_list.end(); }
    size_t get_nr_cells() const { return m_ncells; }
    size_t get_nr_cellsx() const { return m_ncellx; }
    size_t get_nr_unique_pairs() const { return m_atom_neighbor_list.size(); }
    size_t get_direct_nr_unique_pairs(const double max_distance, pele::Array<double> x) const;
    size_t get_maximum_nr_unique_pairs(pele::Array<double> x) const;
    void reset(pele::Array<double> coords);
    bool done() const { return m_container_iterator == m_atom_neighbor_list.end(); }
protected:
    void setup();
    void sanity_check();
    void reset_iterator();
    size_t atom2xbegin(const size_t atom_index) const { return m_ndim * atom_index; }
    template <class T> T loop_pow(const T x, int ex) const;
    size_t atom2cell(const size_t i);
    pele::Array<double> cell2coords(const size_t icell) const;
    bool areneighbors(const size_t icell, const size_t jcell) const;
    double get_minimum_corner_distance2(pele::Array<double>& ic, pele::Array<double>& jc) const;
    double get_norm2(pele::Array<double>& x, pele::Array<double>& y) const;
    void build_cell_neighbors_list();
    void build_atom_neighbors_list();
    void build_linked_lists();
};

template<typename distance_policy>
CellIter<distance_policy>::CellIter(pele::Array<double> const coords,
        std::shared_ptr<distance_policy> dist,
        pele::Array<double> const boxv, const double rcut,
        const double ncellx_scale)
    : m_dist(dist),
      m_coords(coords.copy()),
      m_natoms(coords.size() / m_ndim),
      m_rcut(rcut),
      m_initialised(false),
      m_boxv(boxv.copy()),
      m_ncellx(std::max<size_t>(1, static_cast<size_t>(ncellx_scale * m_boxv[0] / rcut))),     //no of cells in one dimension
      m_ncells(std::pow(m_ncellx, m_ndim)),                                                     //total no of cells
      m_rcell(m_boxv[0] / static_cast<double>(m_ncellx)),                                      //size of cell
      m_hoc(m_ncells),                                                                         //head of chain
      m_ll(m_natoms),                                                                          //linked list
      m_xmin(-0.5 * m_boxv[0]),
      m_xmax(0.5 * m_boxv[0])                                     
{
    if (m_boxv.size() != m_ndim) {
        throw std::runtime_error("CellIter::CellIter: distance policy boxv and cell list boxv differ in size");
    }
    if (*std::min_element(m_boxv.data(), m_boxv.data() + m_ndim) < rcut) {
        throw std::runtime_error("CellIter::CellIter: illegal rcut");
    }
    const double boxv_epsilon = 1e-10;
    const double boxv0 = boxv[0];
    for (size_t i = 1; i < boxv.size(); ++i) {
        if (fabs(boxv0 - boxv[i]) > boxv_epsilon) {
            throw std::runtime_error("CellIter::CellIter: illegal input boxv is not for square box");
        }
        if (boxv[i] < 0) {
            throw std::runtime_error("CellIter::CellIter: illegal inout: boxvector");
        }
    }
    if (m_coords.size() != m_ndim * m_natoms) {
        throw std::runtime_error("CellIter::CellIter: illeal coords size");
    }
    if (m_ncellx == 0) {
        throw std::runtime_error("CellIter::CellIter: illegal lattice spacing");
    }
    if (ncellx_scale < 0) {
        throw std::runtime_error("CellIter::CellIter: illegal input");
    }
    setup();
}

template<typename distance_policy>
size_t CellIter<distance_policy>::get_direct_nr_unique_pairs(const double max_distance, pele::Array<double> x) const
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
size_t CellIter<distance_policy>::get_maximum_nr_unique_pairs(pele::Array<double> x) const
{
    const size_t natoms = x.size() / m_ndim;
    return (natoms * (natoms - 1)) / 2;
}

template <typename distance_policy>
void CellIter<distance_policy>::setup()
{
    m_atom_neighbor_list.reserve(m_natoms * (m_natoms - 1) / 2);
    build_cell_neighbors_list();
    reset(m_coords);
    m_initialised = true;
    //sanity_check();
}

template <typename distance_policy>
void CellIter<distance_policy>::sanity_check()
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
        throw std::runtime_error("CellIter::setup: sanity check failed: too few pairs");
    }
    if (nr_unique_pairs_lists > maximum_nr_unique_pairs) {
        throw std::runtime_error("CellIter::setup: sanity check failed: too many pairs");
    }
}

template <typename distance_policy>
void CellIter<distance_policy>::reset_iterator()
{
    m_atom_neighbor_list.clear();
    m_container_iterator = m_atom_neighbor_list.begin();
}

/*re-build linked lists
 * Algorithm 37 page 552 Understanding Molecular Simulation 2nd ed.
 * start by setting head of chain (hoc of size ncells) to -1 (meaning end of chain)
 * then update linked list so that atom i points to the next atom in the chain,
 * obviously this starts from -1 if it is the only element in the chain. If the next
 * atom i is in the same cell, then the hoc for that cell is set to be i
 * and the linked list at position i will point to the index of the previous atom.
 * This is done iteratively for all atoms.
 */
template <typename distance_policy>
void CellIter<distance_policy>::reset(pele::Array<double> coords)
{
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
                throw std::runtime_error("CellIter::reset: coords are incompatible with boxvector");
            }
        }
    }
    reset_iterator();
    build_linked_lists();
    build_atom_neighbors_list();
}

template <typename distance_policy>
template <class T>
T CellIter<distance_policy>::loop_pow(const T x, int ex) const
{
    T result(1);
    while (--ex > -1) {
        result *= x;
    }
    return result;
}

//this function assumes that particles have been already put in box
template <typename distance_policy>
size_t CellIter<distance_policy>::atom2cell(const size_t i)
{
    assert(i < m_natoms);
    size_t icell = 0;
    for(size_t j = 0; j < m_ndim; ++j) {
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

//returns the coordinates to the corner of one of the cells
template <typename distance_policy>
pele::Array<double> CellIter<distance_policy>::cell2coords(const size_t icell) const
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

template <typename distance_policy>
bool CellIter<distance_policy>::areneighbors(const size_t icell, const size_t jcell) const 
{
    if (icell == jcell) {
        return true;
    }
    // Get "lower-left" corners.
    pele::Array<double> icell_coords = cell2coords(icell);
    pele::Array<double> jcell_coords = cell2coords(jcell);
    // Not neccesary, but makes it clearer.
    for (size_t i = 0; i < m_ndim; ++i) {
        icell_coords[i] += 0.5 * m_rcell;
        jcell_coords[i] += 0.5 * m_rcell;
    }
    return get_minimum_corner_distance2(icell_coords, jcell_coords) <= m_rcut * m_rcut;
}

template <typename distance_policy>
double CellIter<distance_policy>::get_minimum_corner_distance2(pele::Array<double>& ic, pele::Array<double>& jc) const
{
    double result = std::numeric_limits<double>::max();
    for (size_t i = 0; i < m_ndim; ++i) {
        pele::Array<double> corner_i = ic.copy();
        corner_i[i] -= 0.5 * m_rcell;
        for (size_t j = 0; j < m_ndim; ++j) {
            pele::Array<double> corner_j = jc.copy();
            corner_j[j] -= 0.5 * m_rcell;
            const double this_distance2 = get_norm2(corner_i, corner_j);
            if (this_distance2 < result) {
                result = this_distance2;
            }
        }
    } 
    return result;
}

template <typename distance_policy>
double CellIter<distance_policy>::get_norm2(pele::Array<double>& x, pele::Array<double>& y) const
{
    double r[m_ndim];
    m_dist->get_rij(r, x.data(), y.data());
    double r2 = 0;
    for (size_t i = 0; i < m_ndim; ++i) {
        r2 += r[i];
    }
    return r2;
}

template <typename distance_policy>
void CellIter<distance_policy>::build_cell_neighbors_list()
{
    for(size_t i = 0; i < m_ncells; ++i) {
        std::vector<size_t> ineighbors;
        for(size_t j = 0; j < m_ncells; ++j) {
            if (areneighbors(i, j)) { //includes istself as a neighbor
                ineighbors.push_back(j);
            }
        }
        assert(ineighbors.size() > 0);
        m_cell_neighbors.push_back(ineighbors);
    }
    m_cell_neighbors.swap(m_cell_neighbors);
    assert(m_cell_neighbors.size() == m_ncells);
}

template <typename distance_policy>
void CellIter<distance_policy>::build_atom_neighbors_list()
{
    for(size_t i = 0; i < m_natoms; ++i) {
        const size_t icell = atom2cell(i);
        assert(icell < m_cell_neighbors.size());
        //loop through all the neighbouring cells of icell
        for (std::vector<size_t>::const_iterator jit = m_cell_neighbors[icell].begin(); jit != m_cell_neighbors[icell].end(); ++jit) {
            const size_t jcell = *jit;
            long int j = m_hoc[jcell];
            while (j > 0) {
                if (j > static_cast<long int>(i)) { //this should avoid double counting (not sure though)
                    std::pair<size_t, size_t> pair(i, j);
                    m_atom_neighbor_list.push_back(pair);
                }
                j = m_ll[j];
            }
        }
    }
}

template <typename distance_policy>
void CellIter<distance_policy>::build_linked_lists()
{
    m_hoc.assign(-1); //set head of chains to -1 (empty state)
    for(size_t i = 0; i < m_natoms; ++i) {
        size_t icell = atom2cell(i);
        m_ll[i] = m_hoc[icell];
        m_hoc[icell] = i;
    }
}

} // namespace pele

#endif // #ifndef PELE_NEIGHBOR_ITERATOR_H
