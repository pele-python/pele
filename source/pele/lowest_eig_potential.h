#ifndef _PELE_LOWEST_EIG_POTENTIAL_H
#define _PELE_LOWEST_EIG_POTENTIAL_H

#include "base_potential.h"
#include "array.h"
#include <algorithm>
#include <vector>

namespace pele {

class Orthogonalize{
public:
    virtual ~Orthogonalize(){};
    virtual void orthogonalize(Array<double>& coords, Array<double>& vector)=0;
};

class OrthogonalizeTranslational : public Orthogonalize{
protected:
    std::vector<pele::Array<double>> _tr_evec;
    size_t _natoms, _bdim, _ndim;
    double _d, _tol;
public:
    OrthogonalizeTranslational(size_t natoms, size_t bdim, double tol=1e-6);
    virtual ~OrthogonalizeTranslational(){}
    virtual inline void orthogonalize(Array<double>& coords, Array<double>& vector);
};

OrthogonalizeTranslational::OrthogonalizeTranslational(size_t natoms, size_t bdim, double tol):
        _natoms(natoms),_bdim(bdim),_ndim(bdim*natoms),_tol(tol)
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
    bool success = true;
    pele::Array<double> dot_prod(_bdim);
    vector /= norm(vector);
    //generally in this loop success will be set to false
    for (size_t i=0; i<_bdim;++i){
      dot_prod[i] = dot(_tr_evec[i],vector);
      if(std::abs(dot_prod[i]) > _tol){success = false;};
    }

    while (success == false)
    {
        success = true;
        for (size_t i=0; i<_bdim;++i)
        {
            for(size_t j=0;j<_ndim;++j)
                vector[j] -= dot_prod[i]*_tr_evec[i][j];
        }
        vector /= norm(vector);
        for (size_t i=0; i<_bdim;++i){
            dot_prod[i] = dot(_tr_evec[i],vector);
            if (std::abs(dot_prod[i]) > _tol){success = false;};
        }
    }
}

/*
 * Lowest Eigenvalue Potential:
 * gd = g(x+d)
 * _v: normal unit vector in arbitrary direction
 * _d: finite difference step
 * */

class LowestEigPotential : public BasePotential {
protected:
    pele::BasePotential * _potential;
    pele::Array<double> _coords, _coordsd, _g, _gd;
    size_t _bdim, _natoms;
    double _d;
    OrthogonalizeTranslational _orthog;
public:
    LowestEigPotential(pele::BasePotential * potential, pele::Array<double> coords, size_t bdim, double d=1e-6);
    virtual ~LowestEigPotential(){}
    virtual double inline get_energy(pele::Array<double> x);
    virtual double inline get_energy_gradient(pele::Array<double> x, pele::Array<double> grad);
    void reset_coords(pele::Array<double> new_coords);
};

/*constructor*/

LowestEigPotential::LowestEigPotential(pele::BasePotential * potential, pele::Array<double> coords, size_t bdim, double d):
        _potential(potential), _coords(coords), _coordsd(coords.size()), _g(_coords.size()), _gd(_coords.size()), _bdim(bdim),
        _natoms(_coords.size()/_bdim), _d(d), _orthog(_natoms,_bdim)
    {
        double e = _potential->get_energy_gradient(_coords,_g);
    }

/* calculate energy from distance squared, r0 is the hard core distance, r is the distance between the centres */
double inline LowestEigPotential::get_energy(pele::Array<double> x) {

    _orthog.orthogonalize(_coords, x); //takes care of orthogonalizing and normalizing x

    for(size_t i=0;i<x.size();++i)
        _coordsd[i] = _coords[i] + _d*x[i];

    double ed = _potential->get_energy_gradient(_coordsd,_gd);
    _gd -= _g;
    double mu = dot(_gd,x)/_d;

    return mu;
}

/* calculate energy and gradient from distance squared, gradient is in g/|rij|, r0 is the hard core distance, r is the distance between the centres */
double inline LowestEigPotential::get_energy_gradient(pele::Array<double> x, pele::Array<double> grad) {

    _orthog.orthogonalize(_coords, x);  //takes care of orthogonalizing and normalizing x

    for(size_t i=0;i<x.size();++i)
        _coordsd[i] = _coords[i] + _d*x[i];

    double ed = _potential->get_energy_gradient(_coordsd,_gd);
    _gd -= _g;
    double mu = dot(_gd,x)/_d;
    for(size_t i=0;i<x.size();++i)
    {
        grad[i] = 2*_gd[i]/_d - 2*mu*x[i];
    }

    return mu;
}

void LowestEigPotential::reset_coords(pele::Array<double> new_coords)
{
    _coords.assign(new_coords);
    double e = _potential->get_energy_gradient(_coords,_g);
}

}

#endif
