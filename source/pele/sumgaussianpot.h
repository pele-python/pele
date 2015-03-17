#ifndef _PELE_SUMGAUSSIANPOT_H
#define _PELE_SUMGAUSSIANPOT_H

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
    virtual double get_energy_gradient(Array<double> x, Array<double>& grad)
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
};

class SumGaussianPot : public BasePotential {
private:
    const size_t m_bdim;
    const Array<double> m_means;
    const Array<double> m_cov_matrix_diags;
    const size_t m_ngauss;
    std::vector<std::shared_ptr<GaussianPot> > m_potentials;
public:
    SumGaussianPot(size_t bdim, Array<double> means, Array<double> cov_matrix_diags)
        : m_bdim(bdim),
          m_means(means.copy()),
          m_cov_matrix_diags(cov_matrix_diags.copy()),
          m_ngauss(means.size() / m_bdim)
    {
        /*std::cout<<"size means.size "<<m_means.size()<<std::endl;
        std::cout<<"size m_ngauss "<<m_ngauss<<std::endl;*/
        if (means.size() != cov_matrix_diags.size()) {
            throw std::runtime_error("SumGaussianPot: illegal input");
        }
        if (means.size() % m_bdim) {
            throw std::runtime_error("SumGaussianPot: illegal input");
        }
        assert(m_means.size() == m_ngauss*m_bdim);
        /*std::cout<<"mean"<<" "<<m_means<<std::endl;
        std::cout<<"cov"<<" "<<m_cov_matrix_diags<<std::endl;*/

        for (size_t i = 0; i < m_ngauss; ++i) {
            Array<double> this_mean(m_bdim);
            Array<double> this_cov_diag(m_bdim);
            for (size_t j = 0; j < m_bdim; ++j) {
                this_mean[j] = m_means[j + i * m_bdim];
                this_cov_diag[j] = m_cov_matrix_diags[j + i * m_bdim];
            }
            /*std::cout<<"this_mean"<<i<<" "<<this_mean<<std::endl;
            std::cout<<"this_cov"<<i<<" "<<this_cov_diag<<std::endl;*/
            m_potentials.push_back(std::make_shared<GaussianPot>(this_mean, this_cov_diag));
        }
        m_potentials.swap(m_potentials);
    }
    virtual double get_energy(Array<double> x)
    {
        double energy = 0;
        size_t i = 0;
        //std::cout<<"size m_potentials"<<m_potentials.size()<<std::endl;
        for (auto& pot : m_potentials){
            //std::cout<<i<<std::endl;
            energy += pot->get_energy(x);
            //++i;
        }
        return energy;
    }
    virtual double get_energy_gradient(Array<double> x, Array<double>& grad)
    {
        double energy = 0;
        grad.assign(0.0);
        for (auto& pot : m_potentials){
            Array<double> tmp(m_bdim);
            energy += pot->get_energy_gradient(x, tmp);
            for (size_t k = 0; k < grad.size(); ++k) {
                grad[k] += tmp[k];
            }
        }
        return energy;
    }
};
    
} // namespace pele

#endif // #ifndef _PELE_SUMGAUSSIANPOT_H
