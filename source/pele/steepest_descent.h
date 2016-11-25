#ifndef _PELE_STEEPEST_DESCENT_H_
#define _PELE_STEEPEST_DESCENT_H_

#include "optimizer.h"

namespace pele {

/**
 * Deterministic steepest descent, for testing purposes only. 
 */
 
class SteepestDescent : public GradientOptimizer {
private:
    const double m_eta_ini;
    double m_eta; // Iterate as x_{n+1} = x_{n} - \eta\grad E(x_{n}) / norm(\grad E(x_{n}))
    double m_eta_min;
    const bool m_verbose;
public:
    SteepestDescent(std::shared_ptr<BasePotential> potential,
        const Array<double>& x0, const double eta=1,
        const double tol=1e-5, const bool verbose=false)
        : GradientOptimizer(potential, x0, tol),
          m_eta_ini(eta),
          m_eta(eta),
          m_eta_min(2 * std::numeric_limits<double>::epsilon()),
          m_verbose(verbose)
    {
        if (m_verbose) {
            std::cout.precision(std::numeric_limits<double>::digits10);
        }
    }
    
    void one_iteration()
    {
        ++iter_number_;
        compute_func_gradient(x_, f_, g_);
        const double energy_before = f_;
        const Array<double> x_old = x_.copy();
        while (true) {
            double normg = norm(g_);
            for (size_t i=0; i<x_.size(); ++i){
                x_[i] -= m_eta*g_[i]/normg;
            }
            update_rms();
            if (f_ <= energy_before || m_eta < 2 * m_eta_min) {
                m_eta = m_eta_ini;
                break;
            }
            x_ = x_old.copy();
            m_eta *= 0.5;
        }
        if (m_verbose) {
            std::cout << iter_number_ << "\t" << f_ << "\t" << rms_ << "\t" << m_eta << "\n";
        }
    }
    
    void update_rms()
    {
        compute_func_gradient(x_, f_, g_);
        rms_ = norm(g_) / std::sqrt(x_.size());
    }
    
    void reset(Array<double>& x0)
    {
        iter_number_ = 0;
        nfev_ = 0;
        m_eta = m_eta_ini;
        x_.assign(x0);
        update_rms();
    }
    
    void set_eta(const double eta)
    {
        m_eta = eta;
    }
    
    void set_eta_min(const double eta_min)
    {
        m_eta_min = eta_min;
    }
};
    
} // namespace pele

#endif // #ifndef _PELE_STEEPEST_DESCENT_H_
