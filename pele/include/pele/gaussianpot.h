#ifndef _PELE_GAUSSIANPOT_H
#define _PELE_GAUSSIANPOT_H

#include <algorithm>

#include "pele/base_potential.h"

namespace pele {
    
class GaussianPot : public BasePotential {
private:
    const Array<double> m_mean;
    const Array<double> m_cov_diag;
    const double m_gauss_prefactor;
    const size_t m_bdim;
    std::vector<double> m_diag_icov;
public:
    GaussianPot(Array<double> mean, Array<double> cov_diag)
        : m_mean(mean.copy()),
          m_cov_diag(cov_diag.copy()),
          m_gauss_prefactor(-1),
          m_bdim(mean.size())
    {
        assert(m_mean.size() == m_bdim);
        assert(m_cov_diag.size() == m_bdim);
        for (size_t i = 0; i < m_bdim; ++i) {
            m_diag_icov.push_back(1. / m_cov_diag[i]);
        }
        m_diag_icov.shrink_to_fit();
    }
    virtual double get_energy(Array<double> x)
    {
        double xTAx = 0;
        for (size_t i = 0; i < m_bdim; ++i) {
            const double tmp = x[i] - m_mean[i];
            xTAx += tmp * m_diag_icov[i] * tmp;
        }
        return m_gauss_prefactor * std::exp(-0.5 * xTAx);
    }
    virtual double get_energy_gradient(Array<double> x, Array<double> grad)
    {
        double xTAx = 0;
        for (size_t i = 0; i < m_bdim; ++i) {
            const double tmp = x[i] - m_mean[i];
            xTAx += tmp * m_diag_icov[i] * tmp;
            grad[i] = 2 * m_diag_icov[i] * tmp;
        }
        const double energy = m_gauss_prefactor * std::exp(-0.5 * xTAx);
        grad *= -0.5 * energy; 
        return energy;
    }
    virtual double get_energy_gradient_hessian(Array<double> x, Array<double> grad, Array<double> hess)
    {
        double xTAx = 0;
        for (size_t i = 0; i < m_bdim; ++i) {
            const double tmp = x[i] - m_mean[i];
            xTAx += tmp * m_diag_icov[i] * tmp;
            grad[i] = 2 * m_diag_icov[i] * tmp;
        }
        const double energy = m_gauss_prefactor * std::exp(-0.5 * xTAx);
        grad *= -0.5 * energy;
        for (size_t i = 0; i < m_bdim; ++i) {
            for (size_t j = 0; j < m_bdim; ++j) {
                const double xi = x[i] - m_mean[i];
                const double xj = x[j] - m_mean[j];
                hess[i * m_bdim + j] = (xi * xj * m_diag_icov[j] - (i == j)) * m_diag_icov[i] * energy;
            }
        }
        return energy;
    }
    virtual double add_energy_gradient(Array<double> x, Array<double> grad)
    {
        Array<double> grad_term(grad.size());
        const double energy = get_energy_gradient(x, grad_term);
        grad += grad_term;
        return energy;
    }
    virtual double add_energy_gradient_hessian(Array<double> x, Array<double> grad, Array<double> hess)
    {
        Array<double> grad_term(grad.size());
        Array<double> hess_term(hess.size());
        const double energy = get_energy_gradient_hessian(x, grad_term, hess_term);
        grad += grad_term;
        hess += hess_term;
        return energy;
    }
};
    
} // namespace pele

#endif // #ifndef _PELE_GAUSSIANPOT_H
