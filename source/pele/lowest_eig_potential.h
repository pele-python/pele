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
    pele::Array<double> _coords, _coordsd, _g, _gd;
    size_t _bdim, _nparticles;
    double _d;
    LowestEigPotential(pele::BasePotential * potential, pele::Array<double> coords, size_t bdim, double d=1e-6)
public:
    virtual ~LowestEigPotential(){}
    virtual double inline get_energy(pele::Array<double> x);
    virtual double inline get_energy_gradient(pele::Array<double> x, pele::Array<double> grad);
};

/*constructor*/

LowestEigPotential::LowestEigPotential(pele::BasePotential * potential, pele::Array<double> coords, size_t bdim, double d):
        _potential(potential), _coords(coords), _coordsd(coords), _g(_coords.size()), _gd(_coords.size()), _bdim(bdim), _nparticles(_coords.size()/_bdim), _d(d)
    {
        double e = _potential->get_energy_gradient(_coords,_g);
    }

/* calculate energy from distance squared, r0 is the hard core distance, r is the distance between the centres */
double inline LowestEigPotential::get_energy(pele::Array<double> x) {

    _coordsd += x;
    double ed = _potential->get_energy_gradient(_coordsd,_gd);
    _gd -= _g;
    double mu = dot(_gd,x)/_d;
    _coordsd.assign(coords); //reset coordinates
    return mu;
}

/* calculate energy and gradient from distance squared, gradient is in g/|rij|, r0 is the hard core distance, r is the distance between the centres */
double inline LowestEigPotential::get_energy_gradient(pele::Array<double> x, pele::Array<double> grad) {
    _coordsd += x;
    double ed = _potential->get_energy_gradient(_coordsd,_gd);
    _gd -= _g;
    double mu = dot(_gd,x)/_d;
    for(size_t i=0;i<x.size();++i)
    {
        grad[i] = 2*_gd[i] - 2*mu*x[i];
    }
    _coordsd.assign(coords); //reset coordinates
    return mu;
}

}

#endif
