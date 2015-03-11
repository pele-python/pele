#ifndef _BV_SUMGAUSSIANPOT_H
#define _BV_SUMGAUSSIANPOT_H

#include "pele/base_potential.h"

namespace bv {

class SumGaussianPot : public pele::BasePotential {
public:
    SumGaussianPot(size_t bdim, pele::Array<double> means, pele::Array<double> cov)
        : m_bdim(bdim),
          m_means(means),
          m_cov(cov)
    {
        
    }
    virtual double get_energy(pele::Array<double> x)
    {
        return 0;
    }
    virtual double get_energy_gradient(pele::Array<double> x, pele::Array<double> grad)
    {
        return 0;
    }
private:
    const size_t m_bdim;
    const pele::Array<double> m_means;
    const pele::Array<double> m_cov;
};
    
} // namespace bv

#endif // #ifndef _BV_SUMGAUSSIANPOT_H
