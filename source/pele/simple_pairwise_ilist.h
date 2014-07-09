#ifndef PYGMIN_SIMPLE_PAIRWISE_ILIST_H
#define PYGMIN_SIMPLE_PAIRWISE_ILIST_H

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
 * This class loops though atom pairs, computes the distances and gets the
 * value of the energy and gradient from the class pairwise_interaction.
 * pairwise_interaction is a passed parameter and defines the actual
 * potential function.
 */
template<typename pairwise_interaction, typename distance_policy=cartesian_distance<3> >
class SimplePairwiseNeighborList : public BasePotential
{
protected:
    std::shared_ptr<pairwise_interaction> _interaction;
    std::shared_ptr<distance_policy> _dist;
    std::vector<long int> const _neighbor_list;
    static const size_t _ndim = distance_policy::_ndim;

    SimplePairwiseNeighborList(std::shared_ptr<pairwise_interaction> interaction,
            Array<long int> const & neighbor_list, std::shared_ptr<distance_policy> dist=NULL )
        : _interaction(interaction), 
          _dist(dist),
          _neighbor_list(neighbor_list.begin(), neighbor_list.end())
    {
        if(_dist == NULL) _dist = std::make_shared<distance_policy>();
    }

public:
    virtual ~SimplePairwiseNeighborList() {}

    virtual double get_energy(Array<double> x);
    virtual double get_energy_gradient(Array<double> x, Array<double> grad);
};

template<typename pairwise_interaction, typename distance_policy>
inline double SimplePairwiseNeighborList<pairwise_interaction,
       distance_policy>::get_energy_gradient(Array<double> x, Array<double> grad)
{
    double e=0.;
    double gij, dr[_ndim];
    const size_t nlist = _neighbor_list.size();
    assert(x.size() == grad.size());

    grad.assign(0.);

#ifndef NDEBUG
    const size_t natoms = x.size()/_ndim;
    for (size_t i=0; i<nlist; ++i) {
        assert(_neighbor_list[i] < (long int)natoms);
    }
#endif

    for (size_t i=0; i<nlist; i+=2) {
        size_t atom1 = _neighbor_list[i];
        size_t atom2 = _neighbor_list[i+1];
        size_t i1 = _ndim * atom1;
        size_t i2 = _ndim * atom2;

        for (size_t k=0; k<_ndim; ++k) {
            dr[k] = x[i1+k] - x[i2+k];
        }

        double r2 = 0;
        for (size_t k=0;k<_ndim;++k) {
            r2 += dr[k]*dr[k];
        }

        e += _interaction->energy_gradient(r2, &gij, atom1, atom2);
        for (size_t k=0; k<_ndim; ++k) {
            grad[i1+k] -= gij * dr[k];
        }
        for (size_t k=0; k<_ndim; ++k) {
            grad[i2+k] += gij * dr[k];
        }
    }

    return e;
}

template<typename pairwise_interaction, typename distance_policy>
inline double SimplePairwiseNeighborList<pairwise_interaction,
       distance_policy>::get_energy(Array<double> x)
{
    double e=0.;
    size_t const nlist = _neighbor_list.size();

    for (size_t i=0; i<nlist; i+=2) {
        size_t atom1 = _neighbor_list[i];
        size_t atom2 = _neighbor_list[i+1];
        size_t i1 = _ndim * atom1;
        size_t i2 = _ndim * atom2;
        double dr[_ndim];
        for (size_t k=0; k<_ndim; ++k) {
            dr[k] = x[i1+k] - x[i2+k];
        }
        double r2 = 0;
        for (size_t k=0;k<_ndim;++k) {
            r2 += dr[k]*dr[k];
        }
        e += _interaction->energy(r2, atom1, atom2);
    }

    return e;
}
}

#endif
