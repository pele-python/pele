#ifndef _PELE_ML_COST_ONLINE_H
#define _PELE_ML_COST_ONLINE_H

#include "pele/base_potential_online.h"
#include "pele/meta_pow.h"

namespace pele {
    
class MLCostOnline : public BasePotentialOnline {
    const Array<double> m_data;
public:
    virtual ~MLCostOnline() {}
    MLCostOnline(Array<double> data)
        : BasePotentialOnline(data.size()),
          m_data(data.copy())
    {}
    virtual double get_log_probability(const double datum, const Array<double>& parameters) const = 0;
    double get_energy(Array<double> x, const size_t batch_number)
    {
        if (batch_number >= m_data.size()) {
            throw std::runtime_error("MLCostOnline: illegal batch_number input");
        }
        return -get_log_probability(m_data[batch_number], x);
    }
};

class MLCostOnlineGauss : public MLCostOnline {
public:
    MLCostOnlineGauss(Array<double> data)
        : MLCostOnline(data)
    {}
    double get_log_probability(const double datum, const Array<double>& parameters) const
    {
        return -0.5 * (std::log(2 * M_PI * parameters[1]) + pos_int_pow<2>(datum - parameters[0]) / parameters[1]);
    }
};
    
} // namespace pele

#endif // #ifndef _PELE_ML_COST_ONLINE_H
