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
    const size_t _natoms, _iter;
    const double _rcut;

    NeighborIter(pele::Array<double> coords, double rcut, std::shared_ptr<distance_policy> dist=NULL)
        : _dist(dist),
          _coords(coords.copy()),
          _natoms(coords.size()/_ndim),
          _iter(0),
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
    std::vector<std::vector<double>> _neighbors;

    CellIter(pele::Array<double> coords, double boxl, double rcut, std::shared_ptr<distance_policy> dist=NULL)
        : NeighborIter(coords, rcut, dist),
          _ncellx(floor(boxl/rcut)),            //no of cells in one dimension
          _ncells(std::pow(_ncellx, _ndim)),    //total no of cells
          _boxl(boxl),
          _iboxl(1/boxl),
          _rcell(boxl/_ncellx),                 //size of cell
          _hoc(_ncells),                        //head of chain
          _ll(_natoms),                         //linked list
          _neighbors()                          //empty vector
    {}

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
                x += _boxl;
            }
            icell += floor(x / _rcell) * std::pow(_ncellx,j);
        }

        return icell;
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
                d -= _boxl * round(d/_boxl); // DEBUG: adjust distance for pbc, should be using distance
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
            ineighbors = std::vector<double>();
            for(size_t j=0; j<_ncells;++j){
                if (this->_areneighbors(i,j)){
                    ineighbors.push_back(j);
                }
            }
            _neighbors.push_back(neighbor_list);
        }
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
        double E=0;
        size_t icell = _atom2icell(_iter);

        return _atom_pair;
    }

    virtual bool done()
    {
        return _iter == _natoms;
    }
};


}

#endif
