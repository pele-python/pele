#ifndef _PELE_STOCHASTIC_LEVENBERG_MARQUARDT_H_
#define _PELE_STOCHASTIC_LEVENBERG_MARQUARDT_H_

#include "pele/stochastic_gradient_descent.h"

namespace pele {
    
/**
 * Reference:
 * http://yann.lecun.com/exdb/publis/pdf/lecun-98b.pdf
 */
    
class StochasticDiagonalLevenbergMarquardt : public StochasticGradientDescent {
    Array<double> m_running_diagonal_2nd;
    const double m_epsilon;
    const double m_mu;
    const double m_gamma;
    Array<double> m_g2;
public:
    StochasticDiagonalLevenbergMarquardt(std::shared_ptr<BasePotential> potential,
        const Array<double>& x0, const double epsilon=1,
        const double mu=2, const double gamma=0.1,
        const double tol=1e-5, const size_t seed=42,
        const bool verbose=false)
        : StochasticGradientDescent(potential, x0, 1, tol, seed, verbose),
          m_running_diagonal_2nd(x0.size(), 1),
          m_epsilon(epsilon),
          m_mu(mu),
          m_gamma(gamma),
          m_g2(x0.size(), 1)
    {}
    void one_iteration()
    {
        ++iter_number_;
        compute_func_gradient(x_, f_, g_);
        update_running_diagonal_2nd();
        for (size_t i = 0; i < x_.size(); ++i) {
            x_[i] -= m_epsilon * g_[i] / (m_running_diagonal_2nd[i] + m_mu);
        }
        update_rms();
        if (m_verbose) {
            std::cout << iter_number_ << "\t" << f_ << "\t" << rms_ << "\n";
        }
    }
    void compute_func_gradient(Array<double> x, double& func, Array<double> gradient)
    {
        ++nfev_;
        std::shared_ptr<BasePotentialOnline> pot = std::static_pointer_cast<BasePotentialOnline>(potential_);
        func = pot->get_energy_gradient_gradient2_batch(x, std::uniform_int_distribution<size_t>{0, pot->get_nr_batches() - 1}(m_generator), gradient, m_g2);
    }
    void update_running_diagonal_2nd()
    {
        for (size_t i = 0; i < m_g2.size(); ++i) {
            m_running_diagonal_2nd[i] = (1 - m_gamma) * m_running_diagonal_2nd[i] + m_gamma * m_g2[i];
        }
    }
};
    
} // namespace pele

#endif // #ifndef _PELE_STOCHASTIC_LEVENBERG_MARQUARDT_H_
