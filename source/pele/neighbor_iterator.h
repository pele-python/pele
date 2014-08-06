#ifndef PELE_NEIGHBOR_ITERATOR_H
#define PELE_NEIGHBOR_ITERATOR_H

#include <assert.h>
#include <vector>
#include "base_potential.h"
#include "array.h"
#include "distance.h"
#include <iostream>
#include <memory>

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

template<typename distance_policy=cartesian_distance<3> >
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
        size_t index = 0;

        for(size_t j =0;j<_ndim;++j)
        {
            size_t j1 = _natoms*i + j;
            if (_coords[j1] < 0){
                _coords[j1] += _boxl;
            }
            index += floor(x / _rcell) * std::pow(_ncellx,j);
        }
        return index;
    }

public:
    virtual ~NeighborIter() {}

    /*build linked lists*/
    virtual void reset(pele::Array<double> coords)
    {
        _coords.assign(coords);
        _dist.put_in_box(_coords);
        _hoc.assign(-1); //set head of chains to -1 (empty state)
        for(size_t i=0;i<_natoms;++i)
        {

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
