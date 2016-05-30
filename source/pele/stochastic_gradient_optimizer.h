#ifndef _PELE_STOCHASTIC_GRADIENT_OPTIMIZER_H_
#define _PELE_STOCHASTIC_GRADIENT_OPTIMIZER_H_

#include <random>

#include "base_potential_online.h"
#include "optimizer.h"

namespace pele {
    
class StochasticGradientOptimizer : public GradientOptimizer {
protected:
    std::mt19937_64 m_generator;
public:
    virtual ~StochasticGradientOptimizer() {}
    StochasticGradientOptimizer(std::shared_ptr<BasePotential> potential, const Array<double> x0, double tol=1e-4, const size_t seed=42)
        : GradientOptimizer(potential, x0, tol),
          m_generator(seed)
    {}
    virtual void compute_func_gradient(Array<double> x, double& func, Array<double> gradient)
    {
        ++nfev_;
        std::shared_ptr<BasePotentialOnline> pot = std::static_pointer_cast<BasePotentialOnline>(potential_);
        func = pot->get_energy_gradient_batch(x, std::uniform_int_distribution<size_t>{0, pot->get_nr_batches() - 1}(m_generator), gradient);
    }
};
    
} // namespace pele

#endif // #ifndef _PELE_STOCHASTIC_GRADIENT_OPTIMIZER_H_

