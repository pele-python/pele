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
    size_t _atom_pair[2];
    pele::Array<double> _coords;
    const size_t _natoms;
    const double _rcut;

    NeighborIter(pele::Array<double> coords, double rcut, std::shared_ptr<distance_policy> dist=NULL)
        : _dist(dist),
          _coords(coords.copy()),
          _natoms(coords.size()/_ndim),
          _rcut(rcut)
    {
        if(_dist == NULL) _dist = std::make_shared<distance_policy>();
    }

public:
    virtual ~NeighborIter() {}
    /*return next pair of atoms over which to compute energy interaction*/
    virtual size_t* next() =0;
    virtual bool done() =0;
    virtual void reset(pele::Array<double> coords) =0;
};


/*cell list currently only work with box of equal side lengths*/

template<typename distance_policy=periodic_distance<3> >
class CellIter : public NeighborIter<distance_policy>
{
protected:
    const size_t _natoms, _ncellx, _ncells;
    const double _rcell, _boxl, _iboxl;
    pele::Array<long int> _hoc, _ll;

    CellIter(pele::Array<double> coords, double boxl, double rcut, std::shared_ptr<distance_policy> dist=NULL)
        : NeighborIter(coords, rcut, dist),
          _ncellx(floor(boxl/rcut)),            //no of cells in one dimension
          _ncells(std::pow(_ncellx, _ndim)),    //total no of cells
          _boxl(boxl),
          _iboxl(1/boxl),
          _rcell(boxl/_ncellx),                 //size of cell
          _hoc(_ncells),                         //head of chain
          _ll(_natoms)                           //linked list
    {}

    //return cell index from coordinates
    //this function assumes that particles have been already put in box
    size_t _atom2icell(size_t i)
    {
        size_t icell = 0;

        for(size_t j =0;j<_ndim;++j)
        {
            size_t j1 = _natoms*i + j;
            if (_coords[j1] < 0){
                _coords[j1] += _boxl;
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
     * atom i is the same cell, then the hoc for that cell is set to be the i
     * and the linked list at position i will point to the index of the previous atom.
     * This is done iteratively for all atoms.
     */
    virtual void reset(pele::Array<double> coords)
    {
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
            size_t icell = this->_atom2icell(i);
            _ll[i] = _hoc[icell];
            _hoc[icell] = i;
        }
    }


    /*return next pair of atoms over which to compute energy interaction*/
    size_t* next()
    {

        return _atom_pair;
    }
    virtual bool done();
};


}

#endif
