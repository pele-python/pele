#ifndef PYGMIN_SIMPLE_PAIRWISE_POTENTIAL_H
#define PYGMIN_SIMPLE_PAIRWISE_POTENTIAL_H

#include "base_potential.h"
#include "array.h"
#include "distance.h"
#include <memory>

namespace pele
{
/**
 * Define a base class for potentials with simple pairwise interactions that
 * depend only on magnitude of the atom separation
 *
 * This class loops though atom pairs, computes the distances and get's the
 * value of the energy and gradient from the class pairwise_interaction.
 * pairwise_interaction is a passed parameter and defines the actual
 * potential function.
 */
template<typename pairwise_interaction, 
    typename distance_policy = cartesian_distance<3> >
class SimplePairwisePotential : public BasePotential
{
protected:
    static const size_t _ndim = distance_policy::_ndim;
    std::shared_ptr<pairwise_interaction> _interaction;
    std::shared_ptr<distance_policy> _dist;

    SimplePairwisePotential( std::shared_ptr<pairwise_interaction> interaction,
            std::shared_ptr<distance_policy> dist=NULL) 
        : _interaction(interaction), _dist(dist)
    {
        if(_dist == NULL) _dist = std::make_shared<distance_policy>();
    }

public:
    virtual ~SimplePairwisePotential() 
    {}

    virtual double get_energy(Array<double> x);
    virtual double get_energy_gradient(Array<double> x, Array<double> grad)
    {
        grad.assign(0);
        return add_energy_gradient(x, grad);
    }
    virtual double get_energy_gradient_hessian(Array<double> x, Array<double> grad, Array<double> hess)
    {
        grad.assign(0);
        hess.assign(0);
        return add_energy_gradient_hessian(x, grad, hess);
    }
    virtual double add_energy_gradient(Array<double> x, Array<double> grad);
    virtual double add_energy_gradient_hessian(Array<double> x, Array<double> grad, Array<double> hess);
};

template<typename pairwise_interaction, typename distance_policy>
inline double
SimplePairwisePotential<pairwise_interaction,distance_policy>::add_energy_gradient(
        Array<double> x, Array<double> grad)
{
    const size_t natoms = x.size() / _ndim;
    if (_ndim * natoms != x.size()) {
        throw std::runtime_error("x is not divisible by the number of dimensions");
    }
    if (grad.size() != x.size()) {
        throw std::runtime_error("grad must have the same size as x");
    }

    double e = 0.;
    double gij;
    double dr[_ndim];

//    grad.assign(0.);

    for (size_t atomi=0; atomi<natoms; ++atomi) {
        size_t const i1 = _ndim * atomi;
        for (size_t atomj=0; atomj<atomi; ++atomj) {
            size_t const j1 = _ndim * atomj;

            _dist->get_rij(dr, &x[i1], &x[j1]);

            double r2 = 0;
            for (size_t k=0; k<_ndim; ++k) {
                r2 += dr[k]*dr[k];
            }
            e += _interaction->energy_gradient(r2, &gij, atomi, atomj);
            for (size_t k=0; k<_ndim; ++k) {
                grad[i1+k] -= gij * dr[k];
            }
            for (size_t k=0; k<_ndim; ++k) {
                grad[j1+k] += gij * dr[k];
            }
        }
    }
    return e;
}

template<typename pairwise_interaction, typename distance_policy>
inline double SimplePairwisePotential<pairwise_interaction, distance_policy>::add_energy_gradient_hessian(Array<double> x,
        Array<double> grad, Array<double> hess)
{
    double hij, gij;
    double dr[_ndim];
    const size_t N = x.size();
    const size_t natoms = x.size()/_ndim;
    if (_ndim * natoms != x.size()) {
        throw std::runtime_error("x is not divisible by the number of dimensions");
    }
    if (x.size() != grad.size()) {
        throw std::invalid_argument("the gradient has the wrong size");
    }
    if (hess.size() != x.size() * x.size()) {
        throw std::invalid_argument("the Hessian has the wrong size");
    }
//    hess.assign(0.);
//    grad.assign(0.);

    double e = 0.;
    for (size_t atomi=0; atomi<natoms; ++atomi) {
        int i1 = _ndim*atomi;
        for (size_t atomj=0;atomj<atomi;++atomj){
            int j1 = _ndim*atomj;
            _dist->get_rij(dr, &x[i1], &x[j1]);
            double r2 = 0;
            for (size_t k=0;k<_ndim;++k){r2 += dr[k]*dr[k];}

            e += _interaction->energy_gradient_hessian(r2, &gij, &hij, atomi, atomj);

            for (size_t k=0; k<_ndim; ++k)
                grad[i1+k] -= gij * dr[k];
            for (size_t k=0; k<_ndim; ++k)
                grad[j1+k] += gij * dr[k];

            for (size_t k=0; k<_ndim; ++k){
                //diagonal block - diagonal terms
                double Hii_diag = (hij+gij)*dr[k]*dr[k]/r2 - gij;
                hess[N*(i1+k)+i1+k] += Hii_diag;
                hess[N*(j1+k)+j1+k] += Hii_diag;
                //off diagonal block - diagonal terms
                double Hij_diag = -Hii_diag;
                hess[N*(i1+k)+j1+k] = Hij_diag;
                hess[N*(j1+k)+i1+k] = Hij_diag;
                for (size_t l = k+1; l<_ndim; ++l){
                    //diagonal block - off diagonal terms
                    double Hii_off = (hij+gij)*dr[k]*dr[l]/r2;
                    hess[N*(i1+k)+i1+l] += Hii_off;
                    hess[N*(i1+l)+i1+k] += Hii_off;
                    hess[N*(j1+k)+j1+l] += Hii_off;
                    hess[N*(j1+l)+j1+k] += Hii_off;
                    //off diagonal block - off diagonal terms
                    double Hij_off = -Hii_off;
                    hess[N*(i1+k)+j1+l] = Hij_off;
                    hess[N*(i1+l)+j1+k] = Hij_off;
                    hess[N*(j1+k)+i1+l] = Hij_off;
                    hess[N*(j1+l)+i1+k] = Hij_off;
                }
            }
        }
    }
    return e;
}

template<typename pairwise_interaction, typename distance_policy>
inline double SimplePairwisePotential<pairwise_interaction, distance_policy>::get_energy(Array<double> x)
{
    size_t const natoms = x.size()/_ndim;
    if (_ndim * natoms != x.size()) {
        throw std::runtime_error("x is not divisible by the number of dimensions");
    }
    double e=0.;
    double dr[_ndim];

    for (size_t atomi=0; atomi<natoms; ++atomi) {
        size_t i1 = _ndim*atomi;
        for (size_t atomj=0; atomj<atomi; ++atomj) {
            size_t j1 = _ndim*atomj;
            _dist->get_rij(dr, &x[i1], &x[j1]);
            double r2 = 0;
            for (size_t k=0;k<_ndim;++k) {
                r2 += dr[k]*dr[k];
            }
            e += _interaction->energy(r2, atomi, atomj);
        }
    }
    return e;
}
}

#endif
