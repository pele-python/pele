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
    virtual double get_logp(const double datum, const Array<double>& parameters) const = 0;
    double get_energy(Array<double> x, const size_t index)
    {
        if (index >= m_data.size()) {
            throw std::runtime_error("MLCostOnline: illegal index input");
        }
        return -get_logp(m_data[index], x);
    }
    /*
    double get_energy_ogradient(Array<double> x, const size_t index,
        Array<double> ograd)
    {
        if (x.size() != ograd.size()) {
            throw std::runtime_error("MLCostOnline: illegal input x and grad size");
        }
        ograd.assign(0);
        //////////
        return BasePotentialOnline::get_energy(x);
    }
    */
};

class MLCostOnlineGauss : public MLCostOnline {
public:
    MLCostOnlineGauss(Array<double> data)
        : MLCostOnline(data)
    {}
    double get_logp(const double datum, const Array<double>& parameters) const
    {
        return -0.5 * (std::log(2 * M_PI * parameters[1]) + pos_int_pow<2>(datum - parameters[0]) / parameters[1]);
    }
};
    
} // namespace pele

#endif // #ifndef _PELE_ML_COST_ONLINE_H
