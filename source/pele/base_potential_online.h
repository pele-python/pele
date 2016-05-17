#ifndef PELE_POTENTIAL_ONLINE_H
#define PELE_POTENTIAL_ONLINE_H

#include "pele/base_potential.h"

namespace pele {
    
class BasePotentialOnline : public BasePotential {
protected:
    const size_t m_nr_terms;
public:
    BasePotentialOnline(const size_t nr_terms)
        : m_nr_terms(nr_terms)
    {}
    virtual ~BasePotentialOnline() {}
    virtual double get_energy(Array<double> x)
    {
        double energy = 0;
        for (size_t i = 0; i < m_nr_terms; ++i) {
            energy += get_energy(x, i);
        }
        return energy;
    }
    virtual double get_energy_ogradient(Array<double> x, const size_t index, Array<double> ograd)
    {
        numerical_ogradient(x, index, ograd);
        return get_energy(x);
    }
    virtual double get_energy_ogradient_ogradient2(Array<double>x, const size_t index, Array<double> ograd, Array<double> ograd2)
    {
        const double energy = get_energy_ogradient(x, index, ograd);
        numerical_ogradient2(x, index, ograd2);
        return energy;
    }
    virtual void numerical_ogradient(Array<double> x, const size_t index, Array<double> ograd, const double eps=1e-6)
    {
        if (x.size() != ograd.size()) {
            throw std::invalid_argument("ograd.size() be the same as x.size()");
        }
        if (index >= m_nr_terms) {
            throw std::invalid_argument("illegal index");
        }
        Array<double> xnew(x.copy());
        for (size_t i = 0; i < xnew.size(); ++i) {
            xnew[i] -= eps;
            const double eminus = get_energy(xnew, index);
            xnew[i] += 2 * eps;
            const double eplus = get_energy(xnew, index);
            ograd[i] = (eplus - eminus) / (2 * eps);
            xnew[i] = x[i];
        }
    }
    virtual void numerical_ogradient2(Array<double>x, const size_t index, Array<double> ograd2, const double eps=1e-6)
    {
        if (x.size() != ograd2.size()) {
            throw std::invalid_argument("ograd2.size() be the same as x.size()");
        }
        if (index >= m_nr_terms) {
            throw std::invalid_argument("illegal index");
        }
        const size_t N = x.size();
        Array<double> gplus(N);
        Array<double> gminus(N);
        for (size_t i = 0; i < N; ++i) {
            const double backup = x[i];
            x[i] -= eps;
            get_energy_ogradient(x, index, gminus);
            x[i] += 2 * eps;
            get_energy_ogradient(x, index, gplus);
            x[i] = backup;
        }
        for (size_t i = 0; i < N; ++i) {
            ograd2[i] = (gplus[i] - gminus[i]) / (2 * eps);
        }
    }
    virtual double get_energy(Array<double>x, const size_t index)
    {
        throw std::runtime_error("BasePotentialOnline::get_energy(x, i) must be overloaded");
    }
    size_t get_nr_terms() const
    {
        return m_nr_terms;
    }
};
    
} // namespace pele

#endif // #ifndef PELE_POTENTIAL_ONLINE_H
