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
 * in that case _ncellx should be an array rather than a scalar and the definition
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
    static const size_t _ndim = distance_policy::_ndim;
    std::shared_ptr<distance_policy> _dist;
    pele::Array<double> _coords;
    const size_t _natoms;
    const double _rcut;
    bool _initialised;
    const pele::Array<double> _boxv;
    const size_t _ncellx;
    const size_t _ncells;
    const double _rcell;
    pele::Array<long int> _hoc;
    pele::Array<long int> _ll;
    std::vector<std::vector<size_t> > _cell_neighbors;
    std::vector<std::pair<size_t, size_t> > _atom_neighbor_list;
    const_iterator _container_iterator;
    const double _xmin;
    const double _xmax;
public:
    ~CellIter() {}
    CellIter(pele::Array<double> const coords,
        std::shared_ptr<distance_policy> dist,
        pele::Array<double> const boxv, const double rcut,
        const double ncellx_scale=1.0);
    const_iterator begin() const { return _atom_neighbor_list.begin(); }
    const_iterator end() const { return _atom_neighbor_list.end(); }
    size_t get_nr_cells() const { return _ncells; }
    size_t get_nr_cellsx() const { return _ncellx; }
    size_t get_nr_unique_pairs() const { return _atom_neighbor_list.size(); }
    size_t get_direct_nr_unique_pairs(const double max_distance, pele::Array<double> x) const;
    size_t get_maximum_nr_unique_pairs(pele::Array<double> x) const;
    void reset(pele::Array<double> coords);
    bool done() const { return _container_iterator == _atom_neighbor_list.end(); }
protected:
    void _setup();
    void _sanity_check();
    void _reset_iterator();
    size_t _atom2xbegin(const size_t atom_index) const { return _ndim * atom_index; }
    template <class T> T loop_pow(const T x, int ex) const;
    size_t _atom2cell(const size_t i);
    pele::Array<double> _cell2coords(const size_t icell) const;
    bool _areneighbors(const size_t icell, const size_t jcell) const;
    double _get_minimum_corner_distance2(pele::Array<double>& ic, pele::Array<double>& jc) const;
    double _get_norm2(pele::Array<double>& x, pele::Array<double>& y) const;
    void _build_cell_neighbors_list();
    void _build_atom_neighbors_list();
    void _build_linked_lists();
};

template<typename distance_policy>
CellIter<distance_policy>::CellIter(pele::Array<double> const coords,
        std::shared_ptr<distance_policy> dist,
        pele::Array<double> const boxv, const double rcut,
        const double ncellx_scale)
    : _dist(dist),
      _coords(coords.copy()),
      _natoms(coords.size() / _ndim),
      _rcut(rcut),
      _initialised(false),
      _boxv(boxv.copy()),
      _ncellx(std::max<size_t>(1, static_cast<size_t>(ncellx_scale * _boxv[0] / rcut))),     //no of cells in one dimension
      _ncells(std::pow(_ncellx, _ndim)),                                                     //total no of cells
      _rcell(_boxv[0] / static_cast<double>(_ncellx)),                                      //size of cell
      _hoc(_ncells),                                                                         //head of chain
      _ll(_natoms),                                                                          //linked list
      _xmin(-0.5 * _boxv[0]),
      _xmax(0.5 * _boxv[0])                                     
{
    if (_boxv.size() != _ndim) {
        throw std::runtime_error("CellIter::CellIter: distance policy boxv and cell list boxv differ in size");
    }
    if (*std::min_element(_boxv.data(), _boxv.data() + _ndim) < rcut) {
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
    if (_coords.size() != _ndim * _natoms) {
        throw std::runtime_error("CellIter::CellIter: illeal coords size");
    }
    if (_ncellx == 0) {
        throw std::runtime_error("CellIter::CellIter: illegal lattice spacing");
    }
    if (ncellx_scale < 0) {
        throw std::runtime_error("CellIter::CellIter: illegal input");
    }
    this->_setup();
}

template<typename distance_policy>
size_t CellIter<distance_policy>::get_direct_nr_unique_pairs(const double max_distance, pele::Array<double> x) const
{
    size_t nr_unique_pairs = 0;
    const size_t natoms = x.size() / _ndim;
    for (size_t i = 0; i < natoms; ++i) {
        for (size_t j = i + 1; j < natoms; ++j) {
            double rij[_ndim];
            const double* xi = x.data() + _atom2xbegin(i);
            const double* xj = x.data() + _atom2xbegin(j);
            _dist->get_rij(rij, xi, xj);
            double r2 = 0;
            for (size_t k = 0; k < _ndim; ++k) {
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
    const size_t natoms = x.size() / _ndim;
    return (natoms * (natoms - 1)) / 2;
}

template <typename distance_policy>
void CellIter<distance_policy>::_setup()
{
    _atom_neighbor_list.reserve(_natoms * (_natoms - 1) / 2);
    this->_build_cell_neighbors_list();
    this->reset(_coords);
    _initialised = true;
    //_sanity_check();
}

template <typename distance_policy>
void CellIter<distance_policy>::_sanity_check()
{
    const size_t nr_unique_pairs_lists = get_nr_unique_pairs();
    const size_t nr_unique_pairs_direct = get_direct_nr_unique_pairs(_rcut, _coords);
    const size_t maximum_nr_unique_pairs = get_maximum_nr_unique_pairs(_coords);
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
void CellIter<distance_policy>::_reset_iterator()
{
    _atom_neighbor_list.clear();
    _container_iterator = _atom_neighbor_list.begin();
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
    _coords.assign(coords);
    if (periodic_policy_check<distance_policy>::is_periodic) {
        // distance policy is periodic: put particles "back in box" first
        periodic_distance<_ndim>(_boxv).put_in_box(_coords);
    }
    else {
        // distance policy is not periodic: check that particles are inside box
        for (size_t i = 0; i < _coords.size(); ++i) {
            if (_coords[i] < -0.5 * _boxv[0] || _coords[i] > 0.5 * _boxv[0]) {
                std::cout << "_coords[i]: " << _coords[i] << "\n";
                std::cout << "0.5 * _boxv[0]: " << 0.5 * _boxv[0] << std::endl;
                throw std::runtime_error("CellIter::reset: coords are incompatible with boxvector");
            }
        }
    }
    _reset_iterator();
    _build_linked_lists();
    _build_atom_neighbors_list();
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
size_t CellIter<distance_policy>::_atom2cell(const size_t i)
{
    assert(i < _natoms);
    size_t icell = 0;
    for(size_t j = 0; j < _ndim; ++j) {
        const size_t j1 = _atom2xbegin(i) + j;
        assert(j1 < _coords.size());
        double x = _coords[j1];
        // min is needed in case x == _rcell * _ncellx
        const size_t icell_jpart = std::min<size_t>(_ncellx - 1, static_cast<size_t>(((x - _xmin) / (_xmax - _xmin)) * _ncellx));
        assert(icell_jpart == icell_jpart);
        if (icell_jpart >= _ncellx) {
            std::cout << "x: " << x << std::endl;
            std::cout << "_rcell: " << _rcell << std::endl;
            std::cout << "_ndim: " << _ndim << std::endl;
            std::cout << "_ncellx: " << _ncellx << std::endl;
            std::cout << "icell_jpart: " << icell_jpart << std::endl;
        }
        assert(icell_jpart < _ncellx);
        //icell += icell_jpart * std::pow(_ncellx, j);
        icell += icell_jpart * loop_pow(_ncellx, j);
    }
    assert(icell < _ncells);
    return icell;
}

//returns the coordinates to the corner of one of the cells
template <typename distance_policy>
pele::Array<double> CellIter<distance_policy>::_cell2coords(const size_t icell) const
{
    pele::Array<double> cellcorner(_ndim); //coordinate of cell bottom left corner
    std::vector<double> indexes(_ndim, 0); //this array will store indexes, set to 0
    double index = 0;

    //don't change these loops to size_t or the conditions will not hold
    for(int i = _ndim - 1; i >= 0; --i) {
        index = icell;
        for (int j = _ndim - 1; j >= i; --j) {
            //index -= indexes[j] * std::pow(_ncellx, j);
            index -= indexes[j] * loop_pow(_ncellx, j);
        }
        //indexes[i] = floor(index / std::pow(_ncellx, i));
        indexes[i] = floor(index / loop_pow(_ncellx, i));
        cellcorner[i] = _rcell * indexes[i];
    }

    return cellcorner.copy();
}

template <typename distance_policy>
bool CellIter<distance_policy>::_areneighbors(const size_t icell, const size_t jcell) const 
{
    if (icell == jcell) {
        return true;
    }
    // Get "lower-left" corners.
    pele::Array<double> icell_coords = _cell2coords(icell);
    pele::Array<double> jcell_coords = _cell2coords(jcell);
    // Not neccesary, but makes it clearer.
    for (size_t i = 0; i < _ndim; ++i) {
        icell_coords[i] += 0.5 * _rcell;
        jcell_coords[i] += 0.5 * _rcell;
    }
    return _get_minimum_corner_distance2(icell_coords, jcell_coords) <= _rcut * _rcut;
}

template <typename distance_policy>
double CellIter<distance_policy>::_get_minimum_corner_distance2(pele::Array<double>& ic, pele::Array<double>& jc) const
{
    double result = std::numeric_limits<double>::max();
    for (size_t i = 0; i < _ndim; ++i) {
        pele::Array<double> corner_i = ic.copy();
        corner_i[i] -= 0.5 * _rcell;
        for (size_t j = 0; j < _ndim; ++j) {
            pele::Array<double> corner_j = jc.copy();
            corner_j[j] -= 0.5 * _rcell;
            const double this_distance2 = _get_norm2(corner_i, corner_j);
            if (this_distance2 < result) {
                result = this_distance2;
            }
        }
    } 
    return result;
}

template <typename distance_policy>
double CellIter<distance_policy>::_get_norm2(pele::Array<double>& x, pele::Array<double>& y) const
{
    double r[_ndim];
    _dist->get_rij(r, x.data(), y.data());
    double r2 = 0;
    for (size_t i = 0; i < _ndim; ++i) {
        r2 += r[i];
    }
    return r2;
}

template <typename distance_policy>
void CellIter<distance_policy>::_build_cell_neighbors_list()
{
    for(size_t i = 0; i < _ncells; ++i) {
        std::vector<size_t> ineighbors;
        for(size_t j = 0; j < _ncells; ++j) {
            if (this->_areneighbors(i, j)) { //includes istself as a neighbor
                ineighbors.push_back(j);
            }
        }
        assert(ineighbors.size() > 0);
        _cell_neighbors.push_back(ineighbors);
    }
    _cell_neighbors.swap(_cell_neighbors);
    assert(_cell_neighbors.size() == _ncells);
}

template <typename distance_policy>
void CellIter<distance_policy>::_build_atom_neighbors_list()
{
    for(size_t i = 0; i < _natoms; ++i) {
        const size_t icell = this->_atom2cell(i);
        assert(icell < _cell_neighbors.size());
        //loop through all the neighbouring cells of icell
        for (std::vector<size_t>::const_iterator jit = _cell_neighbors[icell].begin(); jit != _cell_neighbors[icell].end(); ++jit) {
            const size_t jcell = *jit;
            long int j = _hoc[jcell];
            while (j > 0) {
                if (j > static_cast<long int>(i)) { //this should avoid double counting (not sure though)
                    std::pair<size_t, size_t> pair(i, j);
                    _atom_neighbor_list.push_back(pair);
                }
                j = _ll[j];
            }
        }
    }
}

template <typename distance_policy>
void CellIter<distance_policy>::_build_linked_lists()
{
    _hoc.assign(-1); //set head of chains to -1 (empty state)
    for(size_t i = 0; i < _natoms; ++i) {
        size_t icell = this->_atom2cell(i);
        _ll[i] = _hoc[icell];
        _hoc[icell] = i;
    }
}

} // namespace pele

#endif // #ifndef PELE_NEIGHBOR_ITERATOR_H
