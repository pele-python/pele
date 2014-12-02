#ifndef _PELE_HARMONIC_H
#define _PELE_HARMONIC_H

#include <algorithm>
#include <functional>

#include "atomlist_potential.h"
#include "base_potential.h"
#include "distance.h"
#include "simple_pairwise_ilist.h"

namespace pele {

class BaseHarmonic : public BasePotential {
protected:
    virtual void _get_distance(const pele::Array<double>& x)=0;
    pele::Array<double> _origin;
    pele::Array<double> _distance;
    double _k;
    size_t _ndim;
    size_t _nparticles;
    BaseHarmonic(pele::Array<double> origin, double k, size_t ndim)
        : _origin(origin.copy()),
          _distance(origin.size()),
          _k(k),
          _ndim(ndim),
          _nparticles(origin.size() / _ndim)
    {}
public:
    virtual ~BaseHarmonic() {}
    virtual double inline get_energy(pele::Array<double> x);
    virtual double inline get_energy_gradient(pele::Array<double> x, pele::Array<double> grad);
    void set_k(double newk) { _k = newk; }
    double get_k() { return _k; }
};


/* calculate energy from distance squared, r0 is the hard core distance, r is the distance between the centres */
double inline BaseHarmonic::get_energy(pele::Array<double> x) 
{
    double norm2 = 0;
    this->_get_distance(x);
    for(size_t i=0;i<x.size();++i)
        norm2 += _distance[i]*_distance[i];
    return 0.5 * _k * norm2;
}

/* calculate energy and gradient from distance squared, gradient is in g/|rij|, r0 is the hard core distance, r is the distance between the centres */
double inline BaseHarmonic::get_energy_gradient(pele::Array<double> x, pele::Array<double> grad) 
{
    assert(grad.size() == _origin.size());
    double norm2 = 0;
    this->_get_distance(x);
    for(size_t i=0;i<x.size();++i){
        norm2 += _distance[i]*_distance[i];
        grad[i] = _k * _distance[i];
    }
    return 0.5 * _k * norm2;
}

/**
 * Simple Harmonic with cartesian distance
 */
class Harmonic : public BaseHarmonic{
public:
    Harmonic(pele::Array<double> origin, double k, size_t ndim)
        : BaseHarmonic(origin, k, ndim)
    {}

    virtual void inline _get_distance(const pele::Array<double>& x)
    {
        assert(x.size() == _origin.size());
        for(size_t i=0;i<x.size();++i)
        {
            _distance[i] = x[i] - _origin[i];;
        }
    }
};

/**
 * Harmonic with cartesian distance and fixed centre of mass
 */
class HarmonicCOM : public BaseHarmonic{
public:
    HarmonicCOM(pele::Array<double> origin, double k, size_t ndim)
        : BaseHarmonic(origin, k, ndim)
    {}

    virtual void inline _get_distance(const pele::Array<double>& x)
    {
        assert(x.size() == _origin.size());
        pele::Array<double> delta_com(_ndim,0);

        for(size_t i=0;i<_nparticles;++i)
        {
            size_t i1 = i*_ndim;
            for(size_t j=0;j<_ndim;++j){
                double d = x[i1+j] - _origin[i1+j];
                _distance[i1+j] = d;
                delta_com[j] += d;
            }
        }

        delta_com /= _nparticles;

        for(size_t i=0;i<_nparticles;++i)
        {
            size_t i1 = i*_ndim;
            for(size_t j=0;j<_ndim;++j)
                _distance[i1+j] -= delta_com[j];
        }
    }
};


/*
 * redo the harmonic potential in the pair potential framework
 */

struct harmonic_interaction {
    double const m_k;
    harmonic_interaction(double k)
        : m_k(k)
    {}

    /* calculate energy from distance squared */
    double inline energy(double r2, size_t atom_i, size_t atom_j) const
    {
        return 0.5 * m_k * r2;
    }

    /* calculate energy and gradient from distance squared, gradient is in g/|rij| */
    double inline energy_gradient(double r2, double *gij, size_t atom_i, size_t atom_j) const
    {
        *gij = -m_k;
        return 0.5 * m_k * r2;
    }

    double inline energy_gradient_hessian(double r2, double *gij, double *hij,
            size_t atom_i, size_t atom_j) const
    {
        *gij = -m_k;
        *hij = 1.;
        return 0.5 * m_k * r2;
    }
};

/**
 * Pairwise harmonic interaction with loops done using atom lists
 */
class HarmonicAtomList : public AtomListPotential<harmonic_interaction, cartesian_distance<3>> {
public:
    HarmonicAtomList(double k, Array<size_t> atoms1, Array<size_t> atoms2)
        : AtomListPotential<harmonic_interaction, cartesian_distance<3>>(
                std::make_shared<harmonic_interaction>(k),
                std::make_shared<cartesian_distance<3>>(), atoms1, atoms2)
    {}

    HarmonicAtomList(double k, Array<size_t> atoms1)
        : AtomListPotential<harmonic_interaction, cartesian_distance<3>>(
                std::make_shared<harmonic_interaction>(k),
                std::make_shared<cartesian_distance<3>>(), atoms1)
    {}
};

/**
 * Pairwise Harmonic potential with neighbor lists
 */
class HarmonicNeighborList : public SimplePairwiseNeighborList<harmonic_interaction> {
public:
    HarmonicNeighborList(double k, Array<size_t> ilist)
        :  SimplePairwiseNeighborList<harmonic_interaction>(
                std::make_shared<harmonic_interaction>(k), ilist)
    {}
};

} // namespace pele

#endif // #ifndef _PELE_HARMONIC_H
