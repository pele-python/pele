#ifndef _PELE_STEEPEST_DESCENT_H_
#define _PELE_STEEPEST_DESCENT_H_

#include "optimizer.h"

namespace pele {

/**
 * Deterministic steepest descent, for testing purposes only. 
 */
 
class SteepestDescent : public GradientOptimizer {
private:
    double m_eta; // Iterate as x_{n+1} = x_{n} - \eta\grad E(x_{n})
    const bool m_verbose;
public:
    SteepestDescent(std::shared_ptr<BasePotential> potential,
        const Array<double>& x0, const double eta=1e-1,
        const double tol=1e-5, const bool verbose=false)
        : GradientOptimizer(potential, x0, tol),
          m_eta(eta),
          m_verbose(verbose)
    {}
    void one_iteration()
    {
        if (m_verbose) {
            std::cout << "before:\t" << iter_number_ << "\t" << f_ << "\t" << g_ << "\t" << m_eta << "\n";
        }
        ++iter_number_;
        compute_func_gradient(x_, f_, g_);
        const double energy_before = f_;
        const Array<double> x_old = x_;
        while (!stop_criterion_satisfied()) {
            x_ -= m_eta * g_;
            update_rms();
            if (f_ < energy_before) {
                break;
            }
            x_ = x_old;
            if (m_eta < 2 * std::numeric_limits<double>::epsilon()) {
                break;
            }
            m_eta *= 0.5;
            if (m_verbose) {
                std::cout << "reduced eta:\t" << m_eta << "\n";
            }
        }
        if (m_verbose) {
            std::cout << "after:\t" << iter_number_ << "\t" << f_ << "\t" << g_ << "\t" << m_eta << "\n";
        }
    }
    void update_rms()
    {
        compute_func_gradient(x_, f_, g_);
        rms_ = norm(g_) / std::sqrt(x_.size());
    }
    void reset(Array<double>& x0, const double eta=1e-1)
    {
        iter_number_ = 0;
        nfev_ = 0;
        m_eta = eta;
        x_.assign(x0);
        update_rms();
    }
};
    
} // namespace pele

#endif // #ifndef _PELE_STEEPEST_DESCENT_H_
