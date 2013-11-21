#ifndef PYGMIN_SIMPLE_PAIRWISE_ILIST_H
#define PYGMIN_SIMPLE_PAIRWISE_ILIST_H

#include <assert.h>
#include <vector>
#include "base_potential.h"
#include "array.h"
#include "distance.h"
#include <iostream>

using std::cout;

namespace pele
{
    /**
     * This class loops though atom pairs, computes the distances and get's the
     * value of the energy and gradient from the class pairwise_interaction.
     * pairwise_interaction is a passed parameter and defines the actual
     * potential function.
     */
    template<typename pairwise_interaction, typename distance_policy=cartesian_distance >
    class SimplePairwiseNeighborList : public BasePotential
    {
        protected:
            pairwise_interaction *_interaction;
            distance_policy *_dist;
            std::vector<long int> const _neighbor_list;

            SimplePairwiseNeighborList(pairwise_interaction *interaction, 
                    Array<long int> const & neighbor_list, distance_policy *dist=NULL ) : 
                _interaction(interaction), 
                _dist(dist),
                _neighbor_list(neighbor_list.begin(), neighbor_list.end())
            {
                if(_dist == 0) _dist = new distance_policy;
            }

        public:
            virtual ~SimplePairwiseNeighborList()
            { 
                if (_interaction != NULL) delete _interaction; 
                if (_dist != NULL) delete _dist; 
            }

            virtual double get_energy(Array<double> x);
            virtual double get_energy_gradient(Array<double> x, Array<double> grad);
    };

    template<typename pairwise_interaction, typename distance_policy>
    inline double SimplePairwiseNeighborList<pairwise_interaction,
           distance_policy>::get_energy_gradient(Array<double> x, Array<double>
                   grad)
    {
        double e=0.;
        double gij, dr[3];
        const size_t n = x.size();
        const size_t natoms = x.size()/3;
        const size_t nlist = _neighbor_list.size();
        assert(x.size() == grad.size());

        for(size_t i=0; i<n; ++i)
            grad[i] = 0.;

#ifndef NDEBUG
        for(size_t i=0; i<nlist; ++i)
            assert(_neighbor_list[i] < (long int)natoms);
#endif

        for(size_t i=0; i<nlist; i+=2) {
            size_t i1 = 3*_neighbor_list[i];
            size_t i2 = 3*_neighbor_list[i+1];

            for(size_t k=0; k<3; ++k)
                dr[k] = x[i1+k] - x[i2+k];

            double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];

            e += _interaction->energy_gradient(r2, &gij);
            for(size_t k=0; k<3; ++k)
                grad[i1+k] -= gij * dr[k];
            for(size_t k=0; k<3; ++k)
                grad[i2+k] += gij * dr[k];
        }

        return e;
    }

    template<typename pairwise_interaction, typename distance_policy>
    inline double SimplePairwiseNeighborList<pairwise_interaction,
           distance_policy>::get_energy(Array<double> x)
    {
        double e=0.;
//        size_t const natoms = x.size()/3;
        size_t const nlist = _neighbor_list.size();

        for(size_t i=0; i<nlist; i+=2) {
            size_t i1 = 3*_neighbor_list[i];
            size_t i2 = 3*_neighbor_list[i+1];
            double dr[3];
            for(size_t k=0; k<3; ++k)
                dr[k] = x[i1+k] - x[i2+k];
            double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
            e += _interaction->energy(r2);
        }

        return e;
    }
}

#endif
