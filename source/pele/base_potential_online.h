#ifndef PELE_POTENTIAL_ONLINE_H
#define PELE_POTENTIAL_ONLINE_H

#include "pele/base_potential.h"

namespace pele {
    
class BasePotentialOnline : public BasePotential {
protected:
    const size_t m_nr_batches;
public:
    BasePotentialOnline(const size_t nr_batches)
        : m_nr_batches(nr_batches)
    {}
    virtual ~BasePotentialOnline() {}
    virtual double get_energy(Array<double> x)
    {
        double energy = 0;
        for (size_t i = 0; i < m_nr_batches; ++i) {
            energy += get_energy(x, i);
        }
        return energy;
    }
    virtual double get_energy_gradient_batch(Array<double> x, const size_t batch_number, Array<double> grad)
    {
        numerical_gradient_batch(x, batch_number, grad);
        return get_energy(x);
    }
    virtual double get_energy_gradient_gradient2_batch(Array<double>x, const size_t batch_number, Array<double> grad,
                                                       Array<double> grad2)
    {
        const double energy = get_energy_gradient_batch(x, batch_number, grad);
        numerical_gradient2_batch(x, batch_number, grad2);
        return energy;
    }
    virtual void numerical_gradient_batch(Array<double> x, const size_t batch_number, Array<double> grad, const double eps=1e-6)
    {
        if (x.size() != grad.size()) {
            throw std::invalid_argument("grad.size() must be the same as x.size()");
        }
        if (batch_number >= m_nr_batches) {
            throw std::invalid_argument("illegal batch_number");
        }
        Array<double> xnew(x.copy());
        for (size_t i = 0; i < xnew.size(); ++i) {
            xnew[i] -= eps;
            const double eminus = get_energy(xnew, batch_number);
            xnew[i] += 2 * eps;
            const double eplus = get_energy(xnew, batch_number);
            grad[i] = (eplus - eminus) / (2 * eps);
            xnew[i] = x[i];
        }
    }
    virtual void numerical_gradient2_batch(Array<double>x, const size_t batch_number, Array<double> grad2, const double eps=1e-6)
    {
        if (x.size() != grad2.size()) {
            throw std::invalid_argument("grad2.size() be the same as x.size()");
        }
        if (batch_number >= m_nr_batches) {
            throw std::invalid_argument("illegal batch_number");
        }
        const size_t N = x.size();
        Array<double> gplus(N);
        Array<double> gminus(N);
        for (size_t i = 0; i < N; ++i) {
            const double backup = x[i];
            x[i] -= eps;
            get_energy_gradient_batch(x, batch_number, gminus);
            x[i] += 2 * eps;
            get_energy_gradient_batch(x, batch_number, gplus);
            x[i] = backup;
        }
        for (size_t i = 0; i < N; ++i) {
            grad2[i] = (gplus[i] - gminus[i]) / (2 * eps);
        }
    }
    virtual double get_energy(Array<double>x, const size_t batch_number)
    {
        throw std::runtime_error("BasePotentialOnline::get_energy(x, i) must be overloaded");
    }
    size_t get_nr_batches() const
    {
        return m_nr_batches;
    }
};
    
} // namespace pele

#endif // #ifndef PELE_POTENTIAL_ONLINE_H
