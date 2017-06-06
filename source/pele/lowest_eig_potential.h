#ifndef _PELE_LOWEST_EIG_POTENTIAL_H
#define _PELE_LOWEST_EIG_POTENTIAL_H

#include "base_potential.h"
#include "array.h"
#include <algorithm>
#include <vector>
#include <cmath>
#include <memory>

namespace pele {

inline void zero_modes_translational(std::vector<pele::Array<double> > & zev,
        size_t natoms, size_t bdim)
{
    double v = 1 / sqrt(natoms);
    size_t N = natoms * bdim;
    for(size_t i=0; i<bdim; ++i) {
        pele::Array<double> evec(N, 0); //initialize array of zeros
        for(size_t j=i; j<N; j+=bdim) {
            evec[j] = v;
        }
        zev.push_back(evec);
    }
}

class Orthogonalize{
public:
    virtual ~Orthogonalize(){};
    virtual void orthogonalize(Array<double> const & coords, Array<double>& vector) =0;
};

class OrthogonalizeTranslational : public Orthogonalize{
protected:
    std::vector<pele::Array<double>> _tr_evec;
    size_t _natoms, _bdim, _ndim;
    double _d, _tol;
public:

    //OrthogonalizeTranslational(size_t natoms, size_t bdim, double tol=1e-6);
    OrthogonalizeTranslational(size_t natoms, size_t bdim, double tol=1e-6)
        : _natoms(natoms), _bdim(bdim), _ndim(bdim*natoms), _tol(tol)
    {
        /*initialize translational eigenvectors to canonical orthonormal basis*/
        zero_modes_translational(_tr_evec, _natoms, _bdim);
    }

    virtual ~OrthogonalizeTranslational() {}

    //virtual inline void orthogonalize(Array<double> const & coords, Array<double>& vector);
    virtual inline void orthogonalize(Array<double> const & coords, Array<double>& vector)
    {
        bool success = true;
        pele::Array<double> dot_prod(_bdim);
        vector /= norm(vector);
        //generally in this loop success will be set to false
        for (size_t i=0; i<_bdim;++i) {
          dot_prod[i] = dot(_tr_evec[i],vector);
          if(std::abs(dot_prod[i]) > _tol){success = false;};
        }

        while (success == false) {
            success = true;
            for (size_t i=0; i<_bdim;++i) {
                for(size_t j=0;j<_ndim;++j) {
                    vector[j] -= dot_prod[i]*_tr_evec[i][j];
                }
            }
            vector /= norm(vector);
            for (size_t i=0; i<_bdim;++i) {
                dot_prod[i] = dot(_tr_evec[i],vector);
                if (std::abs(dot_prod[i]) > _tol) {
                    success = false;
                };
            }
        }
    }

}; //class OrthogonalizeTranslational

/*
 * Lowest Eigenvalue Potential:
 * gd = g(x+d)
 * _v: normal unit vector in arbitrary direction
 * _d: finite difference step
 * */

class LowestEigPotential : public BasePotential {
protected:
    std::shared_ptr<pele::BasePotential> _potential;
    pele::Array<double> _coords, _coordsd, _g, _gd;
    size_t _bdim, _natoms;
    double _d;
    OrthogonalizeTranslational _orthog;
    pele::Array<double> x_opt; //!< Points to the coordinates used in the optimizer
public:

    //LowestEigPotential(std::shared_ptr<pele::BasePotential> potential, pele::Array<double>
    //        coords, size_t bdim, double d=1e-6);
    /*constructor*/
    LowestEigPotential(std::shared_ptr<pele::BasePotential> potential, pele::Array<double>coords,
            size_t bdim, double d=1e-6)
        : _potential(potential), _coords(coords.copy()), _coordsd(coords.size()),
          _g(_coords.size()), _gd(_coords.size()), _bdim(bdim),
          _natoms(_coords.size()/_bdim), _d(d), _orthog(_natoms,_bdim),
          x_opt(coords)
    {
        _potential->get_energy_gradient(_coords,_g);
    }

    void set_x_opt(pele::Array<double> x) {
        x_opt = pele::Array<double>(x);
    }


    virtual ~LowestEigPotential(){}


    //virtual double inline get_energy(pele::Array<double> const & x);
    /* calculate energy from distance squared, r0 is the hard core distance, r is the distance between the centres */
    virtual double inline get_energy(pele::Array<double> const & x)
    {
        if (x_opt != x) {
            throw std::runtime_error("LowestEigPotential::get_energy: x_opt is different from x. Use set_x_opt to pass a pointer.");
        }
        _orthog.orthogonalize(_coords, x_opt); //takes care of orthogonalizing and normalizing x_opt

        for (size_t i = 0; i < x_opt.size(); i++) {
            _coordsd[i] = _coords[i] + _d * x_opt[i];
        }

        _potential->get_energy_gradient(_coordsd, _gd);
        _gd -= _g;
        double mu = dot(_gd,x_opt) / _d;

        return mu;
    }

    //virtual double inline get_energy_gradient(pele::Array<double> x,
    //        pele::Array<double> grad);
    /* calculate energy and gradient from distance squared, gradient is in g/|rij|, r0 is the hard core distance, r is the distance between the centres */
    virtual double inline get_energy_gradient(pele::Array<double> const & x, pele::Array<double> & grad)
    {
        if (x_opt != x) {
            throw std::runtime_error("LowestEigPotential::get_energy_gradient: x_opt is different from x. Use set_x_opt to pass a pointer.");
        }
        _orthog.orthogonalize(_coords, x_opt);  //takes care of orthogonalizing and normalizing x_opt

        for (size_t i = 0; i < x_opt.size(); i++) {
            _coordsd[i] = _coords[i] + _d * x_opt[i];
        }

        _potential->get_energy_gradient(_coordsd, _gd);
        _gd -= _g;
        double mu = dot(_gd, x_opt) / _d;
        for (size_t i = 0; i < x_opt.size(); i++) {
            grad[i] = 2 * _gd[i] / _d - 2 * mu * x[i];
        }

        return mu;
    }

    //void reset_coords(pele::Array<double> new_coords);
    void reset_coords(pele::Array<double> const & new_coords)
    {
        _coords.assign(new_coords);
        _potential->get_energy_gradient(_coords,_g);
    }


};

}//namespace pele

#endif//#ifndef _PELE_LOWEST_EIG_POTENTIAL_H
