#ifndef PELE_POTENTIAL_ONLINE_H
#define PELE_POTENTIAL_ONLINE_H

namespace pele {
    
class BasePotentialOnline {
protected:
    const size_t m_max_index;
public:
    BasePotentialOnline(const size_t max_index)
        : m_max_index(max_index)
    {}
    virtual ~BasePotentialOnline() {}
    virtual double get_energy(Array<double> x)
    {
        double energy = 0;
        for (size_t i = 0; i < m_max_index; ++i) {
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
        // TODO: Add numerical og2
        // TODO: Add remaining functionality (and sgd) and test on simple examples
    }
    virtual double get_energy(Array<double>x, const size_t index)
    {
        throw std::runtime_error("BasePotentialOnline::get_energy(x, i) must be overloaded");
    }
};
    
} // namespace pele

#endif // #ifndef PELE_POTENTIAL_ONLINE_H
