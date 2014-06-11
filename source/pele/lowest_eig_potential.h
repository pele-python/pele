#ifndef _PELE_LOWESTEIGPOTENTIAL_H
#define _PELE_LOWESTEIGPOTENTIAL_H

#include "base_potential.h"
#include "distance.h"
#include <algorithm>
#include <functional>

namespace pele {

/*
 * Lowest Eigenvalue Potential:
 * gd = g(x+d)
 * _v: normal unit vector in arbitrary direction
 * _d: finite difference step
 * */

class LowestEigPotential : public BasePotential {
protected:
    pele::BasePotential * _potential;
    pele::Orthogonalize * _gramschmidt;
    pele::Array<double> _coords, _coordsd, _g, _gd;
    size_t _bdim, _natoms;
    double _d;
public:
    LowestEigPotential(pele::BasePotential * potential, pele::Array<double> coords, size_t bdim, double d=1e-6);
    virtual ~LowestEigPotential(){}
    virtual double inline get_energy(pele::Array<double> x);
    virtual double inline get_energy_gradient(pele::Array<double> x, pele::Array<double> grad);
};

/*constructor*/

LowestEigPotential::LowestEigPotential(pele::BasePotential * potential, pele::Array<double> coords, size_t bdim, double d):
        _potential(potential), _coords(coords), _coordsd(coords), _g(_coords.size()), _gd(_coords.size()), _bdim(bdim),
        _natoms(_coords.size()/_bdim), _d(d)
    {
        double e = _potential->get_energy_gradient(_coords,_g);
    }

/* calculate energy from distance squared, r0 is the hard core distance, r is the distance between the centres */
double inline LowestEigPotential::get_energy(pele::Array<double> x) {

    double normx = norm(x);
    for(size_t i=0;i<x.size();++i)
        _coordsd[i] = _coords[i] + _d*x[i]/normx;

    double ed = _potential->get_energy_gradient(_coordsd,_gd);
    _gd -= _g;
    double mu = dot(_gd,x)/_d;

    return mu;
}

/* calculate energy and gradient from distance squared, gradient is in g/|rij|, r0 is the hard core distance, r is the distance between the centres */
double inline LowestEigPotential::get_energy_gradient(pele::Array<double> x, pele::Array<double> grad) {

    double normx = norm(x);
    for(size_t i=0;i<x.size();++i)
        _coordsd[i] = _coords[i] + _d*x[i]/normx;;

    double ed = _potential->get_energy_gradient(_coordsd,_gd);
    _gd -= _g;
    double mu = dot(_gd,x)/_d;
    for(size_t i=0;i<x.size();++i)
    {
        grad[i] = 2*_gd[i]/_d - 2*mu*x[i];
    }

    return mu;
}


class Orthogonalize{
public:
    virtual ~OrthogonalizeTranslational(){}
    virtual void orthogonalize(Array<double>& coords, Array<double>& vector) = 0;
};


class OrthogonalizeTranslational : public Orthogonalize{
protected:
    std::list<pele::Array<double>> _tr_evec;
    size_t _natoms, _bdim, _ndim;
    double _d;
public:
    OrthogonalizeTranslational(size_t natoms, size_t bdim);
    virtual ~OrthogonalizeTranslational(){}
    virtual inline void orthogonalize(Array<double>& coords, Array<double>& vector);
};

OrthogonalizeTranslational::OrthogonalizeTranslational(size_t natoms, size_t bdim):
        _natoms(natoms),_bdim(bdim),_ndim(bdim*natoms)
{
    /*initialize translational eigenvectors to canonical orthonormal basis*/
    double v = 1/sqrt(_natoms);
    for(size_t i=0;i<_bdim;++i)
    {
        pele::Array<double> evec(_ndim,0); //initialize array of zeros
        for(size_t j=i;j<_ndim;j+=_bdim)
            evec[j] = v;
        _tr_evec.push_back(evec.copy());
    }
}

inline void OrthogonalizeTranslational::orthogonalize(Array<double>& coords, Array<double>& vector)
{
    pele::Array<double> dot_prod(_bdim);
    for(size_t i=0;i<_bdim;++i)
    {
        dot_prod[i] = dot(_tr_evec[i],vector);
    }

    for(size_t i=0;i<_bdim;++i)
    {
        for(size_t j=0;j<_ndim;++j)
        {
            vector[j] -= dot_prod[i]*_tr_evec[i][j];
        }
    }
}




}

#endif
