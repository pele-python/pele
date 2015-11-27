#ifndef _P_SPIN_SPHERICAL_H_
#define _P_SPIN_SPHERICAL_H_

#include "array.h"
#include "base_potential.h"

#include <iostream>
#include <cassert>
#include "pele/meta_pow.h"
#include "pele/combination.h"

namespace pele {

template <size_t p>
class MeanFieldPSpinSpherical : public BasePotential {
protected:
    pele::Array<double> m_interactions, m_spins;
    pele::Array<size_t> m_indexes;
    size_t m_N;
    static const size_t m_p = p;
    pele::combination_generator<size_t*> m_combination_generator;
    pele::combination_generator<size_t*> m_combination_generator_grad;
    size_t m_get_index(size_t* comb){
        size_t idx = comb[0];
        size_t d = m_N;
        for (size_t i = 1; i < m_p; ++i){
            idx += comb[i] * d;
            d *= m_N;
        }
        return idx;
    }
public:
    virtual ~MeanFieldPSpinSpherical() {}
    //right now interactions_ is of size N**p which is not obtimal has it contains all the permutations of the interactions
    //this is fast but costs a lot in terms of memory
    MeanFieldPSpinSpherical(pele::Array<double> interactions_, size_t nspins_)
            : m_interactions(interactions_.copy()),
              m_spins(nspins_),
              m_indexes(nspins_),
              m_N(nspins_),
              m_combination_generator(m_indexes.begin(), m_indexes.end(), m_p),
              m_combination_generator_grad(m_indexes.begin(), m_indexes.end(), m_p-1)
        {
            static_assert(m_p >= 2, "MeanFieldPSpinSpherical, p cannot be less than 2");
            assert((size_t) std::pow(m_N,m_p) == m_interactions.size()); //interactions should actually be of size (N r)
            for (size_t i=0; i<m_N; ++i){
                m_indexes[i] = i;
            }
        }
    virtual inline double get_energy(pele::Array<double> x);
    virtual inline double add_energy_gradient(pele::Array<double> x, pele::Array<double> grad);
    virtual inline double add_energy_gradient_hessian(pele::Array<double> x, pele::Array<double> grad, pele::Array<double> hess);
};

//x has to be of size m_N+1, to include the lagrange multiplier
template <size_t p>
inline double MeanFieldPSpinSpherical<p>::get_energy(pele::Array<double> x){
        assert(x.size() == m_N+1);

        double e = 0;
        size_t comb[m_p];
        bool go_on = true;
        while(go_on){
            go_on = m_combination_generator(comb);
            double sigmaprod = 1;
            for(size_t i=0; i<m_p; ++i){
                sigmaprod *= x[comb[i]];
            }
            e -= m_interactions[this->m_get_index(comb)] * sigmaprod;
        }
        std::copy(x.data(),x.data()+m_N, m_spins.data());
        e += x[m_N] * (dot(m_spins, m_spins) - m_N);
        return e;
}

template <size_t p>
inline double MeanFieldPSpinSpherical<p>::add_energy_gradient(pele::Array<double> x, pele::Array<double> grad){
    size_t comb_full[m_p];
    size_t combg[m_p-1];

    for (size_t i=0; i<m_N; ++i){
        double g = 0;
        bool go_on = true;
        while(go_on){
            go_on = m_combination_generator(combg);
            bool repeated_idx = std::find(combg, combg+m_p-1, i) != combg+m_p-1;
            if (!repeated_idx){
                std::copy(combg, combg+m_p-1, comb_full);
                comb_full[m_p-1] = i;
                double sigmaprod = 1;
                for(size_t j=0; j<m_p-1; ++j){
                    sigmaprod *= x[combg[j]];
                }
                g -= m_interactions[this->m_get_index(comb_full)] * sigmaprod;
            }
        }
        //now set the gradient element i
        grad[i] = g + 2 * x[m_N] * x[i];
    }

    return this->get_energy(x);
}

template <size_t p>
inline double MeanFieldPSpinSpherical<p>::add_energy_gradient_hessian(Array<double> x, Array<double> grad, Array<double> hess){
    size_t comb_full[m_p];
    size_t combh[m_p-2];

    if (m_p == 2){
        hess.assign(0);
    }
    else{
        /*this handles m_p = 2 correctly*/
        for (size_t i=0; i<m_N; ++i){
            for (size_t j=i+1; j<m_N; ++j){
                double h = 0;
                bool go_on = true;
                while(go_on){
                    go_on = m_combination_generator(combh);
                    bool repeated_idx = (std::find(combh, combh+m_p-2, i) != combh+m_p-2 ||
                            std::find(combh, combh+m_p-2, j) != combh+m_p-2);
                    if (!repeated_idx){
                        std::copy(combh, combh+m_p-2, comb_full);
                        comb_full[m_p-2] = i;
                        comb_full[m_p-1] = j;
                        double sigmaprod = 1;
                        for(size_t k=0; k<m_p-2; ++k){
                            sigmaprod *= x[combh[k]];
                        }
                        h -= m_interactions[this->m_get_index(comb_full)] * sigmaprod;
                    }
                }
                //now set the hessian element ij
                hess[i+m_N*j] = h + 2 * x[m_N] * (double) (i==j);
                hess[j+m_N*i] = hess[i+m_N*j];
            }
        }
    }
    return this->add_energy_gradient(x, grad);
}

}
#endif
