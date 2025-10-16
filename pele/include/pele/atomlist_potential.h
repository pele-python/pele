#ifndef _ATOMLIST_POTENTIAL_H_
#define _ATOMLIST_POTENTIAL_H_

#include "array.h"
#include "base_potential.h"
#include <iostream>

namespace pele {

template<typename pairwise_interaction, typename distance_policy>
class AtomListPotential : public BasePotential
{
protected:
    std::shared_ptr<pairwise_interaction> _interaction;
    std::shared_ptr<distance_policy> _dist;
    Array<size_t> _atoms1;
    Array<size_t> _atoms2;
    bool _one_list;
    static const size_t _ndim = distance_policy::_ndim;

    AtomListPotential(
            std::shared_ptr<pairwise_interaction> interaction,
            std::shared_ptr<distance_policy> dist,
            Array<size_t> & atoms1, Array<size_t> & atoms2) 
        : _interaction(interaction),
          _dist(dist),
          _atoms1(atoms1.copy()),
          _atoms2(atoms2.copy()),
          _one_list(false)
    {}

    AtomListPotential( std::shared_ptr<pairwise_interaction> interaction,
            std::shared_ptr<distance_policy> dist, Array<size_t> & atoms1) 
        : _interaction(interaction),
          _dist(dist),
          _atoms1(atoms1.copy()),
          _atoms2(_atoms1),
          _one_list(true)
    {}


public:

    virtual inline double get_energy(Array<double> x)
    {
        double e=0.;
        size_t jstart = 0;
        double dr[_ndim];

        for(size_t i=0; i<_atoms1.size(); ++i) {
            size_t atom1 = _atoms1[i];
            size_t i1 = _ndim * atom1;
            if (_one_list)
                jstart = i+1;
            for(size_t j=jstart; j<_atoms2.size(); ++j) {
                size_t atom2 = _atoms2[j];
                size_t i2 = _ndim * atom2;

                _dist->get_rij(dr, &x[i1], &x[i2]);
                double r2 = 0;
                for (size_t k=0;k<_ndim;++k){r2 += dr[k]*dr[k];}
                e += _interaction->energy(r2, atom1, atom2);
            }
        }

        return e;

    }

    virtual inline double add_energy_gradient(Array<double> x, Array<double> grad)
    {
        if (x.size() != grad.size()) {
            throw std::invalid_argument("the gradient has the wrong size");
        }
        double e=0.;
        double gij;
        size_t jstart = 0;
        double dr[_ndim];;

        for(size_t i=0; i<_atoms1.size(); ++i) {
            size_t atom1 = _atoms1[i];
            size_t i1 = _ndim * atom1;
            if (_one_list){
                jstart = i+1;
            }
            for(size_t j=jstart; j<_atoms2.size(); ++j) {
                size_t atom2 = _atoms2[j];
                size_t i2 = _ndim * atom2;

                _dist->get_rij(dr, &x[i1], &x[i2]);

                double r2 = 0;
                for (size_t k=0;k<_ndim;++k){r2 += dr[k]*dr[k];}

                e += _interaction->energy_gradient(r2, &gij, atom1, atom2);
                for(size_t k=0; k<_ndim; ++k)
                    grad[i1+k] -= gij * dr[k];
                for(size_t k=0; k<_ndim; ++k)
                    grad[i2+k] += gij * dr[k];
            }
        }

        return e;
    }

    virtual inline double add_energy_gradient_hessian(Array<double> x,
            Array<double> grad, Array<double> hess)
    {
        if (x.size() != grad.size()) {
            throw std::invalid_argument("the gradient has the wrong size");
        }
        if (hess.size() != x.size() * x.size()) {
            throw std::invalid_argument("the Hessian has the wrong size");
        }

        double e=0.;
        double hij, gij;
        size_t jstart = 0;
        double dr[_ndim];
        const size_t N = x.size();

        for(size_t i=0; i<_atoms1.size(); ++i) {
            size_t atom1 = _atoms1[i];
            size_t i1 = _ndim * atom1;

            if (_one_list){
                jstart = i+1;
            }
            for(size_t j=jstart; j<_atoms2.size(); ++j) {
                size_t atom2 = _atoms2[j];
                size_t i2 = _ndim * atom2;

                _dist->get_rij(dr, &x[i1], &x[i2]);
                double r2 = 0;
                for (size_t k=0;k<_ndim;++k){r2 += dr[k]*dr[k];}

                e += _interaction->energy_gradient_hessian(r2, &gij, &hij, atom1, atom2);
                for(size_t k=0; k<_ndim; ++k)
                    grad[i1+k] -= gij * dr[k];
                for(size_t k=0; k<_ndim; ++k)
                    grad[i2+k] += gij * dr[k];


                for (size_t k=0; k<_ndim; ++k){
                    //diagonal block - diagonal terms
                    double Hii_diag = (hij+gij)*dr[k]*dr[k]/r2 - gij;
                    hess[N*(i1+k)+i1+k] += Hii_diag;
                    hess[N*(i2+k)+i2+k] += Hii_diag;
                    //off diagonal block - diagonal terms
                    double Hij_diag = -Hii_diag;
                    hess[N*(i1+k)+i2+k] = Hij_diag;
                    hess[N*(i2+k)+i1+k] = Hij_diag;
                    for (size_t l = k+1; l<_ndim; ++l){
                        //diagonal block - off diagonal terms
                        double Hii_off = (hij+gij)*dr[k]*dr[l]/r2;
                        hess[N*(i1+k)+i1+l] += Hii_off;
                        hess[N*(i1+l)+i1+k] += Hii_off;
                        hess[N*(i2+k)+i2+l] += Hii_off;
                        hess[N*(i2+l)+i2+k] += Hii_off;
                        //off diagonal block - off diagonal terms
                        double Hij_off = -Hii_off;
                        hess[N*(i1+k)+i2+l] = Hij_off;
                        hess[N*(i1+l)+i2+k] = Hij_off;
                        hess[N*(i2+k)+i1+l] = Hij_off;
                        hess[N*(i2+l)+i1+k] = Hij_off;
                    }
                }
            }
        }

        return e;
    }

};
}

#endif
