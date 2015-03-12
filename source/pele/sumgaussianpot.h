#ifndef _PELE_SUMGAUSSIANPOT_H
#define _PELE_SUMGAUSSIANPOT_H

#include "pele/base_potential.h"

namespace pele {
    
class GaussianPot : public BasePotential {
public:
    GaussianPot(Array<double> mean, Array<double> cov_diag)
        : m_mean(mean),
          m_cov_diag(cov_diag),
          m_gauss_prefactor(-1),
          m_bdim(mean.size())
    {
        for (size_t i = 0; i < m_bdim; ++i) {
            m_diag_icov.push_back(1 / m_cov_diag.at(i));
        }
        m_diag_icov.swap(m_diag_icov);
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
        grad = Array<double>(m_bdin, 0);
        for (size_t i = 0; i < m_bdim; ++i) {
            const double tmp = x[i] - m_mean[i];
            xTAx += tmp * m_diag_icov[i] * tmp;
            grad[i] = 2 * m_diag_icov[i] * tmp;
        }
        const double energy = m_gauss_prefactor * std::exp(-0.5 * xTAx);
        grad *= -0.5 * energy; 
        return energy;
    }
private:
    const Array<double> m_mean;
    const Array<double> m_cov_diag;
    const double m_gauss_prefactor;
    Array<double> m_diag_icov;
};

class SumGaussianPot : public BasePotential {
public:
    SumGaussianPot(size_t bdim, Array<double> means, Array<double> cov_matrix_diags)
        : m_bdim(bdim),
          m_means(means),
          m_cov(cov_matrix_diags),
          m_ngauss(means.size() / m_bdim)
    {
        if (means.size() != cov.size()) {
            throw std::runtime_error("SumGaussianPot: illegal input");
        }
        if (means.size() % m_bdim) {
            throw std::runtime_error("SumGaussianPot: illegal input");
        }
        for (size_t i = 0; i < m_ngauss; ++i) {
            Array<double> this_mean(m_bdim);
            Array<double> this_cov_diag(m_bdim);
            for (size_t j = i; j < i + m_bdim; ++j) {
                this_mean[j - i] = m_means[j];
                this_cov_diag[j - i] = m_cov_matrix_diags[j];
            }
            m_potentials.push_back(std::make_shared<GaussianPot>(this_mean, this_cov_diag));
        }
        m_potentials.swap(m_potentials);
    }
    virtual double get_energy(Array<double> x)
    {
        double energy = 0;
        for (std::vector<std::shared_ptr<GaussianPot> >::const_iterator i = m_potentials.begin(); i != m_potentials.end(); ++i) {
            energy += *i->get_energy(x);
        }
        return energy;
    }
    virtual double get_energy_gradient(Array<double> x, Array<double> grad)
    {
        double energy = 0;
        grad = Array<double>(m_bdin, 0);
        for (std::vector<std::shared_ptr<GaussianPot> >::const_iterator i = m_potentials.begin(); i != m_potentials.end(); ++i) {
            Array<double> tmp(m_bdin, 0);
            energy += *i->get_energy_gradient(x, tmp);
            grad += tmp;
        }
        return energy;
    }
private:
    const size_t m_bdim;
    const Array<double> m_means;
    const Array<double> m_cov_matrix_diags;
    const size_t m_ngauss;
    std::vector<std::shared_prt<GaussianPot> > m_potentials;
};
    
} // namespace pele

#endif // #ifndef _PELE_SUMGAUSSIANPOT_H
