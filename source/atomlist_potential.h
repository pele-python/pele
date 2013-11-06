#ifndef _ATOMLIST_POTENTIAL_H_
#define _ATOMLIST_POTENTIAL_H_

#include "array.h"
#include "base_potential.h"
#include <iostream>

namespace pele {
    template<typename pairwise_interaction,
        typename distance_policy>
    class AtomListPotential : public BasePotential
    {
    protected:
        pairwise_interaction * _interaction;
        distance_policy * _dist;
        Array<size_t> _atoms1;
        Array<size_t> _atoms2;
        bool _one_list;

        AtomListPotential(pairwise_interaction * interaction, distance_policy * dist,
                Array<size_t> & atoms1, Array<size_t> & atoms2) :
                    _interaction(interaction),
                    _dist(dist),
                    _atoms1(atoms1.copy()),
                    _atoms2(atoms2.copy()),
                    _one_list(false)
        {
//            std::cout << "using two lists with natoms " << _atoms1.size() << " " << _atoms2.size() << "\n";
        }

        AtomListPotential(pairwise_interaction * interaction, distance_policy * dist,
                Array<size_t> & atoms1) :
                    _interaction(interaction),
                    _dist(dist),
                    _atoms1(atoms1.copy()),
                    _atoms2(_atoms1),
                    _one_list(true)
        {
//            std::cout << "using one list with natoms " << _atoms2.size() << "\n";
        }


    public:
        virtual ~AtomListPotential()
        {
            if (_interaction != NULL) delete _interaction;
            if (_dist != NULL) delete _dist;
            _interaction = NULL;
            _dist = NULL;

        }

        virtual inline double get_energy(Array<double> x)
        {
            double e=0.;
            size_t jstart = 0;

            for(size_t i=0; i<_atoms1.size(); ++i) {
                int i1 = 3*_atoms1[i];
                if (_one_list)
                    jstart = i+1;
                for(size_t j=jstart; j<_atoms2.size(); ++j) {
                    int i2 = 3*_atoms2[j];
                    double dr[3];
                    _dist->get_rij(dr, &x[i1], &x[i2]);
                    double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                    e += _interaction->energy(r2);
                }
            }

            return e;

        }

        virtual inline double add_energy_gradient(Array<double> x, Array<double> grad)
        {
            double e=0.;
            double gij, dr[3];
            size_t jstart = 0;

            for(size_t i=0; i<_atoms1.size(); ++i) {
                int i1 = 3*_atoms1[i];
                if (_one_list){
                    jstart = i+1;
                }
                for(size_t j=jstart; j<_atoms2.size(); ++j) {
                    int i2 = 3*_atoms2[j];

                    _dist->get_rij(dr, &x[i1], &x[i2]);

                    double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                    e += _interaction->energy_gradient(r2, &gij);
                    for(size_t k=0; k<3; ++k)
                        grad[i1+k] -= gij * dr[k];
                    for(size_t k=0; k<3; ++k)
                        grad[i2+k] += gij * dr[k];
                }
            }

            return e;
        }

    };
}

#endif
