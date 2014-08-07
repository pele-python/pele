#ifndef PELE_NEIGHBOR_ITERATOR_H
#define PELE_NEIGHBOR_ITERATOR_H

#include <assert.h>
#include <vector>
#include "base_potential.h"
#include "array.h"
#include "distance.h"
#include <iostream>
#include <memory>
#include <exception>

using std::cout;

namespace pele
{
/**
 * This class returns atom pairs
 */
template<typename distance_policy>
class NeighborIter
{
protected:
    std::shared_ptr<distance_policy> _dist;
    static const size_t _ndim = distance_policy::_ndim;
    size_t _atoms_pair[2];
    pele::Array<double> _coords;
    const size_t _natoms, _iter;
    const double _rcut;
    bool _initialised;

    NeighborIter(pele::Array<double> coords, double rcut, std::shared_ptr<distance_policy> dist=NULL)
        : _dist(dist),
          _coords(coords.copy()),
          _natoms(coords.size()/_ndim),
          _iter(0),
          _rcut(rcut),
          _initialised(false)
    {
        if(_dist == NULL) _dist = std::make_shared<distance_policy>();
    }

    virtual void _setup();
public:
    virtual ~NeighborIter() {}
    /*return next pair of atoms over which to compute energy interaction*/
    virtual size_t* next() =0;
    virtual bool done() =0;
    virtual void reset(pele::Array<double> coords) =0;
};


/*cell list currently only work with box of equal side lengths
 * cell lists are currently not implemented for non cubic boxes:
 * in that case _ncellx should be an array rather than a scalar and the definition
 * of ncells and rcell would change to array. This implies that in all the function these
 * would need to be replace with the correct array element. This adds room for errors
 * so in this first implementation we do not account for that scenario
 * */

template<typename distance_policy=periodic_distance<3> >
class CellIter : public NeighborIter<distance_policy>
{
protected:
    const pele::Array<double> _boxv;
    const size_t _natoms, _ncellx, _ncells, _atomi, _atomj, _neigh, _neigh_max;
    const size_t _neigh; // the count of the neighboring cells
    const double _rcell;
    pele::Array<long int> _hoc, _ll;
    std::vector<std::vector<double>> _neighbors;
    //std::list<std::pair<size_t, size_t> > _atom_neghbor_list;

    CellIter(pele::Array<double> coords, pele::Array<double> boxv, double rcut, std::shared_ptr<distance_policy> dist=NULL)
        : NeighborIter(coords, rcut, dist),
          _boxv(boxv),
          _ncellx(floor(boxv[0]/rcut)),            //no of cells in one dimension
          _ncells(std::pow(_ncellx, _ndim)),    //total no of cells
          _atomi(-1),
          _atomj(-1),
          _neigh(0),
          _neigh_max(0),
          _rcell(boxv[0]/_ncellx),                 //size of cell
          _hoc(_ncells),                        //head of chain
          _ll(_natoms),                         //linked list
          _neighbors()                          //empty vector
    {
        try{
            pele::Array<double> dp_boxv(distance_policy::_box);
            if(dp_boxv != _boxv){
                throw std::runtime_error("distance policy boxv and cell list boxv differ in size");
            }
        }
        catch (exception& e){
            cout << "this type of distance does not have a box array. Exception: " << e.what() << '\n';
        }
        if(_boxv[0] != _boxv[_ndim-1] || _boxv[0] != _boxv[_ndim]){
            throw std::runtime_error("cell lists not implemented for non cubic box");
        }
        this->_setup();
    }

    //returns the coordinates to the corner of one of the cells
    double * _cell2coords(size_t icell)
    {
        double cellcorner[_ndim]; //coordinate of cell bottom left corner
        std::vector<double> indexes(_ndim,0); //this array will store indexes, set to 0
        double index = 0;

        for(size_t i = _ndim - 1; i >= 0; --i)
        {
            index = icell;
            for (int j = _ndim - 1; j >= i; --j)
            {
                index -= indexes[j] * std::pow(_ncellx,j);
            }

            indexes[i] = floor(index / std::pow(_ncellx,i));
            cellcorner[i] = _rcell * indexes[i];
        }

        return cellcorner;
    }


    //test whether 2 cells are neighbours
    bool _areneighbors(size_t icell, size_t jcell)
    {
        double icell_coords[_ndim];
        double jcell_coords[_ndim];

        icell_coords = this->_cell2coords(icell);
        jcell_coords = this->_cell2coords(jcell);
        //compute difference

        for (int i=0;i<_ndim;++i){
            double dxmin;
            bool dxmin_trial = false;
            icell_coords[i] -= jcell_coords[i];

            for(int j=0;j<=1;++j){ //DEBUG should include j=-1 like in jake's implementation?
                double d = icell_coords[i] + j*_rcell;
                d -= _boxv[0] * round(d/_boxv[0]); // DEBUG: adjust distance for pbc, assuming regular cubic box
                if (std::abs(d) < dxmin || !dxmin_trial){
                    dxmin = d;
                    dxmin_trial = true;
                }
            }
            icell_coords[i] =  dxmin;
        }

        double r2 = 0;
        for (int i=0;i<_ndim;++i){
            r2 +=  icell_coords[i]*icell_coords[i];
        }

        return r2 <= _rcut*_rcut;
    }

    void _build_cell_neighbors_list()
    {
        for(size_t i=0; i<_ncells;++i){
            std::vector<double> ineighbors;
            for(size_t j=0; j<_ncells;++j){
                if (this->_areneighbors(i,j)){
                    ineighbors.push_back(j);
                }
            }
            _neighbors.push_back(ineighbors);
        }
    }

    void _setup()
    {
        this->_build_cell_neighbors_list();
        this->reset(_coords);
        _initialised = true;
    }

    //return cell index from coordinates
    //this function assumes that particles have been already put in box
    size_t _atom2cell(size_t i)
    {
        size_t icell = 0;

        for(size_t j =0;j<_ndim;++j)
        {
            size_t j1 = _natoms*i + j;
            double x = _coords[j1];

            if (x < 0){
                x += _boxv[j];
            }
            icell += floor(x / _rcell) * std::pow(_ncellx,j);
        }

        return icell;
    }


public:
    virtual ~NeighborIter() {}

    /*re-build linked lists
     * Algorithm 37 page 552 Understanding Molecular Simulation 2nd ed.
     * start by setting head of chain (hoc of size ncells) to -1 (meaning end of chain)
     * then update linked list so that atom i points to the next atom in the chain,
     * obviously this starts from -1 if it is the only element in the chain. If the next
     * atom i is in the same cell, then the hoc for that cell is set to be i
     * and the linked list at position i will point to the index of the previous atom.
     * This is done iteratively for all atoms.
     */
    virtual void reset(pele::Array<double> coords)
    {
        _atomi=-1;
        _atomj=-1;
        _coords.assign(coords);

        try{
            _dist.put_in_box(_coords);
        }
        catch (exception& e){
            cout << "put in box not implemented for this type of distance. Exception: " << e.what() << '\n';
        }

        _hoc.assign(-1); //set head of chains to -1 (empty state)

        for(size_t i=0;i<_natoms;++i)
        {
            size_t icell = this->_atom2cell(i);
            _ll[i] = _hoc[icell];
            _hoc[icell] = i;
        }
    }

    /*return next pair of atoms over which to compute energy interaction*/
    size_t* next()
    {
        size_t icell; // the cell that atomi is in
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
                icell = this->_atom2cell(_atomi);
                _neigh_max = _neighbors[icell].size();
                _neigh = 0;
            } else{
                // we're done with jcell
                jcell = _neighbors[icell][_neigh];
                ++_neigh;
            }
            _atomj = _hoc(jcell);
        }

        _atoms_pair[0] = _atomi;
        _atoms_pair[1] = _atomj;
        ++_iter;

        return _atom_pair;
    }

    virtual bool done() const
    {
        return _iter == _natoms*_natoms;
    }
};


}

#endif
