#ifndef _PELE_STOCHASTIC_GRADIENT_DESCENT_H_
#define _PELE_STOCHASTIC_GRADIENT_DESCENT_H_

#include "stochastic_gradient_optimizer.h"

namespace pele {
    
/**
 * Stochastic gradient descent with fixed scalar stepsize.
 */
 
class StochasticGradientDescent : public StochasticGradientOptimizer {
    const double m_eta;
protected:
    const bool m_verbose;
public:
    virtual ~StochasticGradientDescent() {}
    
    StochasticGradientDescent(std::shared_ptr<BasePotential> potential,
        const Array<double>& x0, const double eta=1,
        const double tol=1e-5, const size_t seed=42,
        const bool verbose=false)
        : StochasticGradientOptimizer(potential, x0, tol, seed),
          m_eta(eta),
          m_verbose(verbose)
    {
        if (m_verbose) {
            std::cout.precision(std::numeric_limits<double>::digits10);
        }
    }
    
    virtual void one_iteration()
    {
        ++iter_number_;
        compute_func_gradient(x_, f_, g_);
        const double se = std::sqrt(m_eta);
        x_ /= se;
        x_ -= se * g_;
        x_ *= se;
        update_rms();
        if (m_verbose) {
            std::cout << iter_number_ << "\t" << f_ << "\t" << rms_ << "\n";
        }
    }
    
    void update_rms()
    {
        rms_ = norm(g_) / std::sqrt(x_.size());
    }
    
    void reset(Array<double>& x0)
    {
        iter_number_ = 0;
        nfev_ = 0;
        x_.assign(x0);
        update_rms();
    }
};
    
} // namespace pele

#endif // #ifndef _PELE_STOCHASTIC_GRADIENT_DESCENT_H_

