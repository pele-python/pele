#ifndef _PELE_STOCHASTIC_GRADIENT_OPTIMIZER_H_
#define _PELE_STOCHASTIC_GRADIENT_OPTIMIZER_H_

#include <random>

#include "base_potential_online.h"
#include "optimizer.h"

namespace pele {
    
class StochasticGradientOptimizer : public GradientOptimizer {
protected:
    std::shared_ptr<BasePotentialOnline> m_online_potential;
    std::mt19937_64 m_generator;
    std::uniform_int_distribution<size_t> m_rand;
public:
    virtual ~StochasticGradientOptimizer() {}
    StochasticGradientOptimizer(std::shared_ptr<BasePotential> potential, const Array<double> x0,
                                double tol=1e-4, const size_t seed=42)
        : GradientOptimizer(potential, x0, tol),
          m_online_potential(std::static_pointer_cast<BasePotentialOnline>(potential_)),
          m_generator(seed),
          m_rand(0, m_online_potential->get_nr_batches()-1)
    {}

    virtual void compute_func_gradient(Array<double> x, double& func, Array<double> gradient)
    {
        ++nfev_;
        const size_t batch_number = m_rand(m_generator);
        func = m_online_potential->get_energy_gradient_batch(x, batch_number, gradient);
    }
};
    
} // namespace pele

#endif // #ifndef _PELE_STOCHASTIC_GRADIENT_OPTIMIZER_H_

