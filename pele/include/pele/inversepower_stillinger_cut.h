#ifndef _PELE_INVERSEPOWER_STILLIGNER_CUT_H
#define _PELE_INVERSEPOWER_STILLIGNER_CUT_H

#include "pele/simple_pairwise_potential.h"
#include "pele/cell_list_potential.h"

namespace pele {
    
struct InversePowerStillinger_cut_interaction {
    // Inverse power law potential, see JCP 83, 4767 (1984),
    // http://dx.doi.org/10.1063/1.449840
    const double m_pow, m_rcut, m_rcut2, m_q0, m_q1, m_q2;
    const pele::Array<double> m_radii;
    InversePowerStillinger_cut_interaction(const size_t pow, const pele::Array<double> radii, const double rcut)
        : m_pow(pow),
          m_rcut(rcut),
          m_rcut2(m_rcut*m_rcut),
          m_q0(-0.5*(m_pow+1)*(m_pow+2)/std::pow(m_rcut, m_pow)),
          m_q1(m_pow*(m_pow+2)/std::pow(m_rcut, m_pow+1)),
          m_q2(-0.5*m_pow*(m_pow+1)/std::pow(m_rcut, m_pow+2)),
          m_radii(radii.copy())
    {
        if (radii.size() == 0) {
            throw std::runtime_error("InversePowerStillingerCut: illegal input: radii");
        }
    }
    double energy(double r2, size_t atomi, size_t atomj) const
    {
        if (r2 > m_rcut2){
            return 0.;
        }
        const double a = m_radii[atomi] + m_radii[atomj];
        const double an = power_inp2(a*a, m_pow);
        const double irn = 1./power_inp2(r2, m_pow);
        const double E = an * (irn + m_q0 + m_q1*std::sqrt(r2) + m_q2*r2);
        return E;
    }
    // calculate energy and gradient from distance squared, gradient is in -(dv/drij)/|rij|
    double energy_gradient(double r2, double* gij, size_t atomi, size_t atomj) const
    {
        if (r2 > m_rcut2){
            *gij = 0;
            return 0.;
        }
        const double a = m_radii[atomi] + m_radii[atomj];
        const double an = power_inp2(a*a, m_pow);
        const double r = std::sqrt(r2);
        const double irn = 1./this->power(r, m_pow);
        const double E = an * (irn + m_q0 + m_q1*r + m_q2*r2);
        *gij = an * (m_pow*irn/r2 - m_q1/r - 2*m_q2);
        return E;
    }
    double energy_gradient_hessian(double r2, double* gij, double* hij, size_t atomi, size_t atomj) const
    {
        if (r2 > m_rcut2){
            *gij = 0;
            *hij = 0;
            return 0.;
        }
        const double a = m_radii[atomi] + m_radii[atomj];
        const double an = power_inp2(a*a, m_pow);
        const double r = std::sqrt(r2);
        const double irn = 1./this->power(r, m_pow);
        const double E = an * (irn + m_q0 + m_q1*r + m_q2*r2);
        *gij = an * (m_pow*irn/r2 - m_q1/r - 2*m_q2);
        *hij = an * (m_pow*(m_pow+1)*irn/r2 + 2*m_q2);
        return E;
    }

    // Compute inp ** n.
    // See Skiena, p. 48.
    template<class T>
    T power(const T inp, const size_t n) const
    {
        if (n == 0) {
            return 1;
        }
        const T x = power(inp, n / 2);
        if (n % 2) {
            return inp * x * x;
        }
        else {
            return x * x;
        }
    }
    // Compute inp ** n, given inp ** 2.
    template<class T>
    T power_inp2(const T inp2, const size_t n) const
    {
        const T x = power(inp2, n / 2);
        if (n % 2) {
            return std::sqrt(inp2) * x;
        }
        else {
            return x;
        }
    }
};

template<size_t ndim>
class InversePowerStillingerCut : public SimplePairwisePotential<InversePowerStillinger_cut_interaction, cartesian_distance<ndim> > {
public:
    InversePowerStillingerCut(const size_t pow, const pele::Array<double> radii, double rcut)
        : SimplePairwisePotential<InversePowerStillinger_cut_interaction, cartesian_distance<ndim> >(
            std::make_shared<InversePowerStillinger_cut_interaction>(pow, radii, rcut),
            std::make_shared<cartesian_distance<ndim> >())
    {}
};

template<size_t ndim>
class InversePowerStillingerCutPeriodic : public SimplePairwisePotential<InversePowerStillinger_cut_interaction, periodic_distance<ndim> > {
public:
    InversePowerStillingerCutPeriodic(const size_t pow, const pele::Array<double> radii, double rcut, const pele::Array<double> boxvec)
        : SimplePairwisePotential<InversePowerStillinger_cut_interaction, periodic_distance<ndim> >(
            std::make_shared<InversePowerStillinger_cut_interaction>(pow, radii, rcut),
            std::make_shared<periodic_distance<ndim> >(boxvec))
    {}
};

template<size_t ndim>
class InversePowerStillingerCutPeriodicCellLists : public CellListPotential<InversePowerStillinger_cut_interaction, periodic_distance<ndim> > {
public:
    InversePowerStillingerCutPeriodicCellLists(const size_t pow, const pele::Array<double> radii, double rcut,
    const pele::Array<double> boxvec, double ncellx_scale)
        : CellListPotential<InversePowerStillinger_cut_interaction, periodic_distance<ndim> >(
            std::make_shared<InversePowerStillinger_cut_interaction>(pow, radii, rcut),
            std::make_shared<periodic_distance<ndim> >(boxvec),
            boxvec, rcut, ncellx_scale)
    {}
};
    
} // namespace pele

#endif // #ifndef _PELE_INVERSEPOWER_STILLIGNER_CUT_H
