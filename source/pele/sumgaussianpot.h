#ifndef _PELE_SUMGAUSSIANPOT_H
#define _PELE_SUMGAUSSIANPOT_H

#include "pele/base_potential.h"

namespace pele {

class SumGaussianPot : public BasePotential {
public:
    SumGaussianPot(size_t bdim, Array<double> means, Array<double> cov)
        : m_bdim(bdim),
          m_means(means),
          m_cov(cov)
    {
        
    }
    virtual double get_energy(Array<double> x)
    {
        return 0;
    }
    virtual double get_energy_gradient(Array<double> x, Array<double> grad)
    {
        return 0;
    }
private:
    const size_t m_bdim;
    const Array<double> m_means;
    const Array<double> m_cov;
};
    
} // namespace pele

#endif // #ifndef _PELE_SUMGAUSSIANPOT_H
