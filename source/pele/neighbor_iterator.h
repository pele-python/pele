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

namespace pele
{

/*cell list currently only work with box of equal side lengths
 * cell lists are currently not implemented for non cubic boxes:
 * in that case _ncellx should be an array rather than a scalar and the definition
 * of ncells and rcell would change to array. This implies that in all the function these
 * would need to be replace with the correct array element. This adds room for errors
 * so in this first implementation we do not account for that scenario
 * */

template<typename distance_policy = periodic_distance<3> >
class CellIter
{
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
    pele::Array<long int> _hoc, _ll;
    std::vector<std::vector<double> > _cell_neighbors;
    std::vector<std::pair<size_t, size_t> > _atom_neighbor_list;
    const_iterator _container_iterator;

public:

    CellIter(pele::Array<double> coords, pele::Array<double> boxv, const double rcut, const double ncellx_scale = 1.0)
        : _dist(std::make_shared<distance_policy>(boxv.copy())), //DEBUG: this won't work in principle with all distance_policies
          _coords(coords.copy()),
          _natoms(coords.size() / _ndim),
          _rcut(rcut),
          _initialised(false),
          _boxv(boxv.copy()),
          _ncellx(ncellx_scale * floor(boxv[0] / rcut)),            //no of cells in one dimension
          _ncells(std::pow(_ncellx, _ndim)),      //total no of cells
          _rcell(boxv[0] / _ncellx),                //size of cell
          _hoc(_ncells),                          //head of chain
          _ll(_natoms),                          //linked list
          _cell_neighbors(),                      //empty vector
          _atom_neighbor_list()
    {
        if (_boxv.size() != _ndim) {
            throw std::runtime_error("CellIter::CellIter: distance policy boxv and cell list boxv differ in size");
        }
        if (*std::min_element(_boxv.data(), _boxv.data() + _ndim) < rcut) {
            throw std::runtime_error("CellIter::CellIter: illegal rcut");
        }
        this->_setup();
    }

    ~CellIter() {}

    const_iterator begin() { return _atom_neighbor_list.begin(); }
    const_iterator end() { return _atom_neighbor_list.end(); }
    size_t get_nr_unique_pairs() const { return _atom_neighbor_list.size(); }

    void _setup()
    {
        this->_build_cell_neighbors_list();
        this->reset(_coords);
        _initialised = true;
    }

    void _reset_iterator()
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

    void reset(pele::Array<double> coords)
    {
        _coords.assign(coords);
        _dist->put_in_box(_coords);

        this->_reset_iterator();
        this->_build_linked_lists();
        this->_build_atom_neighbors_list();
    }

    bool done() const
    {
        return _container_iterator == _atom_neighbor_list.end();
    }

    //return cell index from coordinates
    //this function assumes that particles have been already put in box
    inline size_t _atom2cell(const size_t i)
    {
        size_t icell = 0;
        for(size_t j = 0; j < _ndim; ++j) {
            size_t j1 = _natoms * i + j;
            double x = _coords[j1];
            if (x < 0) {
                x += _boxv[j];
            }
            icell += floor(x / _rcell) * std::pow(_ncellx, j);
        }
        return icell;
    }

    //returns the coordinates to the corner of one of the cells
    pele::Array<double> _cell2coords(const size_t icell)
    {
        pele::Array<double> cellcorner(_ndim); //coordinate of cell bottom left corner
        std::vector<double> indexes(_ndim, 0); //this array will store indexes, set to 0
        double index = 0;

        //don't change these loops to size_t or the conditions will not hold
        for(int i = _ndim - 1; i >= 0; --i) {
            index = icell;
            for (int j = _ndim - 1; j >= i; --j) {
                index -= indexes[j] * std::pow(_ncellx, j);
            }

            indexes[i] = floor(index / std::pow(_ncellx, i));
            cellcorner[i] = _rcell * indexes[i];
        }

        return cellcorner;
    }


    //test whether 2 cells are neighbours
    bool _areneighbors(size_t icell, size_t jcell)
    {
        pele::Array<double> icell_coords(this->_cell2coords(icell));
        pele::Array<double> jcell_coords(this->_cell2coords(jcell));
        //compute difference

        for (size_t i = 0; i < _ndim; ++i) {
            double dxmin;
            bool dxmin_trial = false;
            icell_coords[i] -= jcell_coords[i];

            for(size_t j = 0; j <= 1; ++j) { //DEBUG should include j=-1 like in jake's implementation?
                double d = icell_coords[i] + j * _rcell;
                d -= _boxv[0] * round(d / _boxv[0]); // DEBUG: adjust distance for pbc, assuming regular cubic box
                if (std::abs(d) < dxmin || !dxmin_trial) {
                    dxmin = d;
                    dxmin_trial = true;
                }
            }
            icell_coords[i] =  dxmin;
        }

        double r2 = 0;
        for (size_t i = 0; i < _ndim; ++i) {
            r2 +=  icell_coords[i] * icell_coords[i];
        }

        return r2 <= _rcut * _rcut;
    }

    void _build_cell_neighbors_list()
    {
        for(size_t i = 0; i < _ncells; ++i) {
            std::vector<double> ineighbors;
            for(size_t j = 0; j < _ncells; ++j) {
                if (this->_areneighbors(i, j)) { //includes istself as a neighbor
                    ineighbors.push_back(j);
                }
            }
            _cell_neighbors.push_back(ineighbors);
        }
    }

    inline void _build_atom_neighbors_list()
    {
        for(size_t i = 0; i < _natoms; ++i) {
             const size_t icell = this->_atom2cell(i);
             //loop through all the neighbouring cells of icell
             for(auto& jcell : _cell_neighbors[icell]) {
                 double j = _hoc[jcell];
                 while (j > 0) {
                     if (j > i) { //this should avoid double counting (not sure though)
                         std::pair<size_t, size_t> pair(i, j);
                         _atom_neighbor_list.push_back(pair);
                     }
                     j = _ll[j];
                 }
             }
        }
    }

    inline void _build_linked_lists()
    {
        _hoc.assign(-1); //set head of chains to -1 (empty state)

        for(size_t i = 0; i < _natoms; ++i) {
            size_t icell = this->_atom2cell(i);
            _ll[i] = _hoc[icell];
            _hoc[icell] = i;
        }
    }


    /*const size_t _neigh, _neigh_max; // the count of the neighboring cells
          * _atomi(-1),
           _atomj(-1),
    */
    /*return next pair of atoms over which to compute energy interaction
    size_t* next()
    {
        // _icell; the cell that atomi is in [this needs to be a member]
        size_t jcell; // the cell thta atomj is in

        if (_iter >= _natoms*_natoms){
            throw std::runtime_error("_iter exceeds 2*_natoms");
        }

        if (_atomi > 0){
            _atomj = _ll[_atomj];
        }
        // if atomi < 0 then atomj will also be less than zero (it's the first iteration)

        while (_atomj < 0){
            // we're at the end of a linked list (or in the first iteration)
            if (_neigh >= _neigh_max){
                // we're done with atomi, we need to move to the next atom.
                ++_atomi;
                _icell = this->_atom2cell(_atomi);
                _neigh_max = _cell_neighbors[_icell].size();
                _neigh = 0;
            } else{
                // we're done with jcell
                jcell = _cell_neighbors[_icell][_neigh];
                ++_neigh;
            }
            _atomj = _hoc(jcell);
        }

        _atoms_pair[0] = _atomi;
        _atoms_pair[1] = _atomj;
        ++_iter;

        return _atom_pair;
    }*/
};


} //namespace pele

#endif //#ifndef PELE_NEIGHBOR_ITERATOR_H
