#ifndef _PELE_HS_WCA_H
#define _PELE_HS_WCA_H

#include <algorithm>
#include <memory>

#include "atomlist_potential.h"
#include "cell_list_potential.h"
#include "distance.h"
#include "meta_pow.h"
#include "simple_pairwise_ilist.h"
#include "simple_pairwise_potential.h"
#include "base_interaction.h"

namespace pele {

/**
 * This a WCA-like potential, which mostly avoids square roots and
 * extrapolates linearly "into the hard core".
 * That should be useful to minimize HS-WCA-like systems.
 * The pair potential is:
 * V_\text{sfHS-WCA}(r^2) = 0                                               \text{ if } r \geq r_S
 * V_\text{sfHS-WCA}(r^2) = V_\text{fHS-WCA}(r^2)                           \text{ if } r_\times < r < r_S
 * V_\text{sfHS-WCA}(r^2) = E_\times - (\sqrt{r^2} - r_\times)G_\times      \text{ if } r \leq r_\times
 * Here:
 * E_\times = V_\text{fHS-WCA}(r_\times)
 * G_\times = \text{grad}[V_\text{fHS-WCA}](r_\times)
 * And:
 * radius_sum : sum of hard radii
 * r_S = (1 + \alpha) * radius_sum
 * r_\times = radius_sum + \delta
 * The choice of the delta parameter below is somewhat arbitrary and
 * could probably be optimised.
 * Computing the gradient GX at the point where we go from fWCA to
 * linear is somewhat confusing because the gradient is originally
 * computed as grad / (-r).
 */
struct sf_HS_WCA_interaction : BaseInteraction {
    const double m_eps;
    const double m_alpha;
    const double m_delta;
    double m_prfac;
    sf_HS_WCA_interaction(const double eps, const double alpha, const double delta=1e-10)
        : m_eps(eps),
          m_alpha(alpha),
          m_delta(delta)
    {
        double prfac_inner = (2 + m_alpha) * m_alpha;
        m_prfac = (prfac_inner * prfac_inner) * (prfac_inner / std::sqrt(2));
        if (eps < 0) {
            throw std::runtime_error("HS_WCA: illegal input: eps");
        }
        if (alpha < 0) {
            throw std::runtime_error("HS_WCA: illegal input: alpha");
        }
    }

    double energy(const double r2, const double radius_sum) const
    {
        const double r_S = (1 + m_alpha) * radius_sum;
        const double r_S2 = pos_int_pow<2>(r_S);
        if (r2 > r_S2) {
            // Energy: separation larger than soft shell.
            return 0;
        }
        // r2 <= r_S2: we have to compute the remaining quantities.
        const double r_X = radius_sum + m_delta;
        const double r_X2 = pos_int_pow<2>(r_X);
        const double r_H2 = pos_int_pow<2>(radius_sum);
        if (r2 > r_X2) {
            // Energy: separation in fHS-WCA regime.
            const double dr = r2 - r_H2;
            const double ir2 = 1.0 / (dr * dr);
            const double ir6 = ir2 * ir2 * ir2;
            const double C3 = m_prfac * r_H2 * r_H2 * r_H2;
            const double C6 = C3 * C3;
            const double C6_ir6 = C6 * ir6;
            const double C12_ir12 = C6_ir6 * C6_ir6;
            return std::max<double>(0, 4. * m_eps * (-C6_ir6 + C12_ir12) + m_eps);
        }
        // r2 <= r_X2
        // Energy: separation in linear regime.
        const double dr = r_X2 - r_H2;
        const double ir2 = 1.0 / (dr * dr);
        const double ir6 = ir2 * ir2 * ir2;
        const double C3 = m_prfac * r_H2 * r_H2 * r_H2;
        const double C6 = C3 * C3;
        const double C6_ir6 = C6 * ir6;
        const double C12_ir12 = C6_ir6 * C6_ir6;
        const double EX = std::max<double>(0, 4. * m_eps * (-C6_ir6 + C12_ir12) + m_eps);
        const double GX = (m_eps * (- 48. * C6_ir6 + 96. * C12_ir12) / dr) * (-r_X);
        return EX + GX * (std::sqrt(r2) - r_X);
    }

    double energy_gradient(const double r2, double *const gij, const double radius_sum) const
    {
        const double r_S = (1 + m_alpha) * radius_sum;
        const double r_S2 = pos_int_pow<2>(r_S);
        if (r2 > r_S2) {
            // Energy, gradient: separation larger than soft shell.
            *gij = 0;
            return 0;
        }
        // r2 <= r_S2: we have to compute the remaining quantities.
        const double r_X = radius_sum + m_delta;
        const double r_X2 = pos_int_pow<2>(r_X);
        const double r_H2 = pos_int_pow<2>(radius_sum);
        if (r2 > r_X2) {
            // Energy, gradient: separation in fHS-WCA regime.
            const double dr = r2 - r_H2;
            const double ir2 = 1.0 / (dr * dr);
            const double ir6 = ir2 * ir2 * ir2;
            const double C3 = m_prfac * r_H2 * r_H2 * r_H2;
            const double C6 = C3 * C3;
            const double C6_ir6 = C6 * ir6;
            const double C12_ir12 = C6_ir6 * C6_ir6;
            *gij = m_eps * (- 48. * C6_ir6 + 96. * C12_ir12) / dr;
            return std::max<double>(0, 4. * m_eps * (-C6_ir6 + C12_ir12) + m_eps);
        }
        // r2 <= r_X2
        // Energy, gradient: separation in linear regime.
        const double dr = r_X2 - r_H2;
        const double ir2 = 1.0 / (dr * dr);
        const double ir6 = ir2 * ir2 * ir2;
        const double C3 = m_prfac * r_H2 * r_H2 * r_H2;
        const double C6 = C3 * C3;
        const double C6_ir6 = C6 * ir6;
        const double C12_ir12 = C6_ir6 * C6_ir6;
        const double EX = std::max<double>(0, 4. * m_eps * (-C6_ir6 + C12_ir12) + m_eps);
        const double GX = (m_eps * (- 48. * C6_ir6 + 96. * C12_ir12) / dr) * (-r_X);
        *gij = GX / (-r_X);
        return EX + GX * (std::sqrt(r2) - r_X);
    }

    double energy_gradient_hessian(const double r2, double *const gij,
            double *const hij, const double radius_sum) const
    {
        const double r_S = (1 + m_alpha) * radius_sum;
        const double r_S2 = pos_int_pow<2>(r_S);
        if (r2 > r_S2) {
            // Energy, gradient, hessian: separation larger than soft shell.
            *gij = 0;
            *hij = 0;
            return 0;
        }
        // r2 <= r_S2: we have to compute the remaining quantities.
        const double r_X = radius_sum + m_delta;
        const double r_X2 = pos_int_pow<2>(r_X);
        const double r_H2 = pos_int_pow<2>(radius_sum);
        if (r2 > r_X2) {
            // Energy, gradient, hessian: separation in fHS-WCA regime.
            const double dr = r2 - r_H2;
            const double ir2 = 1.0 / (dr * dr);
            const double ir6 = ir2 * ir2 * ir2;
            const double C3 = m_prfac * r_H2 * r_H2 * r_H2;
            const double C6 = C3 * C3;
            const double C6_ir6 = C6 * ir6;
            const double C12_ir12 = C6_ir6 * C6_ir6;
            *gij = m_eps * (- 48. * C6_ir6 + 96. * C12_ir12) / dr;
            *hij = -*gij + m_eps * ( -672. * C6_ir6 + 2496. * C12_ir12)  * r2 * ir2;
            return std::max<double>(0, 4. * m_eps * (-C6_ir6 + C12_ir12) + m_eps);
        }
        // r2 <= r_X2
        // Energy, gradient, hessian: separation in linear regime.
        const double dr = r_X2 - r_H2;
        const double ir2 = 1.0 / (dr * dr);
        const double ir6 = ir2 * ir2 * ir2;
        const double C3 = m_prfac * r_H2 * r_H2 * r_H2;
        const double C6 = C3 * C3;
        const double C6_ir6 = C6 * ir6;
        const double C12_ir12 = C6_ir6 * C6_ir6;
        const double EX = std::max<double>(0, 4. * m_eps * (-C6_ir6 + C12_ir12) + m_eps);
        const double GX = (m_eps * (- 48. * C6_ir6 + 96. * C12_ir12) / dr) * (-r_X);
        *gij = GX / (-r_X);
        *hij = 0;
        return EX + GX * (std::sqrt(r2) - r_X);
    }
    /**
     * This can be used to plot the potential, as evaluated numerically.
     */
    void evaluate_pair_potential(const double rmin, const double rmax,
            const size_t nr_points, const double radius_sum,
            std::vector<double>& x, std::vector<double>& y) const
    {
        x = std::vector<double>(nr_points, 0);
        y = std::vector<double>(nr_points, 0);
        const double rdelta = (rmax - rmin) / (nr_points - 1);
        for (size_t i = 0; i < nr_points; ++i) {
            x.at(i) = rmin + i * rdelta;
            y.at(i) = energy(x.at(i) * x.at(i), radius_sum);
        }
    }
};

/**
 * Fast pairwise interaction for Hard Sphere + Weeks-Chandler-Andersen (fHS_WCA) potential, refer to S. Martiniani CPGS pp 20
 * well depth _eps and scaling factor (shell thickness = sca * R, where R is the hard core radius),
 * sca determines the thickness of the shell
 */
struct HS_WCA_interaction : BaseInteraction {
    const double _eps;
    const double _sca;
    const double _infty;
    const double _prfac;

    HS_WCA_interaction(const double eps, const double sca)
        : _eps(eps),
          _sca(sca),
          _infty(std::pow(10.0, 50)),
          _prfac(0.5 * pos_int_pow<6>((2 + _sca) * _sca))
    {

    }

    /* calculate energy from distance squared, r0 is the hard core distance, r is the distance between the centres */
    double inline energy(const double r2, const double radius_sum) const
    {
        const double r02 = radius_sum * radius_sum;
        if (r2 <= r02) {
            return _infty;
        }
        const double coff = radius_sum * (1.0 + _sca); //distance at which the soft cores are at contact
        if (r2 > coff * coff) {
            return 0;
        }
        const double dr = r2 - r02; // note that dr is the difference of the squares
        const double ir = 1.0 / dr;
        const double r02_ir = r02 * ir;
        const double C_ir_6 = _prfac * pos_int_pow<6>(r02_ir);
        const double C_ir_12 = C_ir_6 * C_ir_6;
        return compute_energy(C_ir_6, C_ir_12);
    }

    /* calculate energy and gradient from distance squared, gradient is in g/|rij|, radius_sum is the hard core distance, r is the distance between the centres */
    double inline energy_gradient(const double r2, double *const gij, const double radius_sum) const
    {
        const double r02 = radius_sum * radius_sum;
        if (r2 <= r02) {
            *gij = _infty;
            return _infty;
        }
        const double coff = radius_sum * (1.0 + _sca); //distance at which the soft cores are at contact
        if (r2 > coff * coff) {
            *gij = 0.;
            return 0.;
        }
        const double dr = r2 - r02; // note that dr is the difference of the squares
        const double ir = 1.0 / dr;
        const double r02_ir = r02 * ir;
        const double C_ir_6 = _prfac * pos_int_pow<6>(r02_ir);
        const double C_ir_12 = C_ir_6 * C_ir_6;
        *gij = _eps * (- 48. * C_ir_6 + 96. * C_ir_12) * ir; //this is -g/r, 1/dr because powers must be 7 and 13
        return compute_energy(C_ir_6, C_ir_12);
    }

    double inline energy_gradient_hessian(const double r2, double *const gij, double *const hij, const double radius_sum) const
    {
        const double r02 = radius_sum * radius_sum;
        if (r2 <= r02) {
            *gij = _infty;
            *hij = _infty;
            return _infty;
        }
        const double coff = radius_sum * (1.0 + _sca); //distance at which the soft cores are at contact
        if (r2 > coff * coff) {
            *gij = 0.;
            *hij = 0.;
            return 0.;
        }
        const double dr = r2 - r02; // note that dr is the difference of the squares
        const double ir = 1.0 / dr;
        const double r02_ir = r02 * ir;
        const double C_ir_6 = _prfac * pos_int_pow<6>(r02_ir);
        const double C_ir_12 = C_ir_6 * C_ir_6;
        *gij = _eps * (- 48. * C_ir_6 + 96. * C_ir_12) * ir; //this is -g/r, 1/dr because powers must be 7 and 13
        *hij = -*gij + _eps * ( -672. * C_ir_6 + 2496. * C_ir_12)  * r2 * ir * ir;
        return compute_energy(C_ir_6, C_ir_12);
    }

    /**
     * If r2 > r02 && r2 <= coff * coff, this computes the energy.
     * Returning the maximum of 0 and the enrgy term makes the result
     * strictly positive also in the regime close to zero where there
     * can be numerical issues.
     */
    double compute_energy(const double C_ir_6, const double C_ir_12) const
    {
        return std::min<double>(_infty, std::max<double>(0, 4. * _eps * (-C_ir_6 + C_ir_12) + _eps));
    }

    /**
     * This can be used to plot the potential, as evaluated numerically.
     */
    void evaluate_pair_potential(const double rmin, const double rmax, const size_t nr_points, const double radius_sum, std::vector<double>& x, std::vector<double>& y) const
    {
        x = std::vector<double>(nr_points, 0);
        y = std::vector<double>(nr_points, 0);
        const double rdelta = (rmax - rmin) / (nr_points - 1);
        for (size_t i = 0; i < nr_points; ++i) {
            x.at(i) = rmin + i * rdelta;
            y.at(i) = energy(x.at(i) * x.at(i), radius_sum);
        }
    }
};

//
// combine the components (interaction, looping method, distance function) into
// defined classes
//

/**
 * Pairwise HS_WCA potential
 */
template<size_t ndim>
class HS_WCA : public SimplePairwisePotential< sf_HS_WCA_interaction, cartesian_distance<ndim> > {
public:
    HS_WCA(double eps, double sca, Array<double> radii)
        : SimplePairwisePotential< sf_HS_WCA_interaction, cartesian_distance<ndim> >(
                std::make_shared<sf_HS_WCA_interaction>(eps, sca),
                radii,
                std::make_shared<cartesian_distance<ndim> >(),
                sca
            )
    {}
};

/**
 * Pairwise HS_WCA potential in a rectangular box
 */
template<size_t ndim>
class HS_WCAPeriodic : public SimplePairwisePotential< sf_HS_WCA_interaction, periodic_distance<ndim> > {
public:
    HS_WCAPeriodic(double eps, double sca, Array<double> radii, Array<double> const boxvec)
        : SimplePairwisePotential< sf_HS_WCA_interaction, periodic_distance<ndim> > (
                std::make_shared<sf_HS_WCA_interaction>(eps, sca),
                radii,
                std::make_shared<periodic_distance<ndim> >(boxvec),
                sca
                )
    {}
};

/**
 * Pairwise HS_WCA potential in a rectangular box with shear
 */
template<size_t ndim>
class HS_WCALeesEdwards : public SimplePairwisePotential< sf_HS_WCA_interaction, leesedwards_distance<ndim> > {
public:
    HS_WCALeesEdwards(double eps, double sca, Array<double> radii, Array<double> const boxvec,
                const double shear)
        : SimplePairwisePotential< sf_HS_WCA_interaction, leesedwards_distance<ndim> > (
                std::make_shared<sf_HS_WCA_interaction>(eps, sca),
                radii,
                std::make_shared<leesedwards_distance<ndim> >(boxvec, shear),
                sca
                )
    {}
};

template<size_t ndim>
class HS_WCACellLists : public CellListPotential< sf_HS_WCA_interaction, cartesian_distance<ndim> > {
public:
    HS_WCACellLists(double eps, double sca, Array<double> radii, Array<double> const boxvec,
            const double ncellx_scale=1.0, const bool balance_omp=true)
    : CellListPotential< sf_HS_WCA_interaction, cartesian_distance<ndim> >(
            std::make_shared<sf_HS_WCA_interaction>(eps, sca),
            std::make_shared<cartesian_distance<ndim> >(),
            boxvec,
            2 * (1 + sca) * *std::max_element(radii.begin(), radii.end()), // rcut
            ncellx_scale, radii, sca, balance_omp)
    {
        if (boxvec.size() != ndim) {
            throw std::runtime_error("HS_WCA: illegal input: boxvec");
        }
    }
    size_t get_nr_unique_pairs() const { return CellListPotential< sf_HS_WCA_interaction, cartesian_distance<ndim> >::m_celliter->get_nr_unique_pairs(); }
};

template<size_t ndim>
class HS_WCAPeriodicCellLists : public CellListPotential< sf_HS_WCA_interaction, periodic_distance<ndim> > {
public:
    HS_WCAPeriodicCellLists(double eps, double sca, Array<double> radii, Array<double> const boxvec,
            const double ncellx_scale=1.0, const bool balance_omp=true)
    : CellListPotential< sf_HS_WCA_interaction, periodic_distance<ndim> >(
            std::make_shared<sf_HS_WCA_interaction>(eps, sca),
            std::make_shared<periodic_distance<ndim> >(boxvec),
            boxvec,
            2 * (1 + sca) * *std::max_element(radii.begin(), radii.end()), // rcut
            ncellx_scale, radii, sca, balance_omp)
    {}
    size_t get_nr_unique_pairs() const { return CellListPotential< sf_HS_WCA_interaction, periodic_distance<ndim> >::m_celliter->get_nr_unique_pairs(); }
};

/**
 * HS_WCA potential using CellLists in a rectangular box with shear
 */
template<size_t ndim>
class HS_WCALeesEdwardsCellLists : public CellListPotential< sf_HS_WCA_interaction, leesedwards_distance<ndim> > {
public:
    HS_WCALeesEdwardsCellLists(double eps, double sca, Array<double> radii, Array<double> const boxvec,
            const double shear, const double ncellx_scale=1.0, const bool balance_omp=true)
    : CellListPotential< sf_HS_WCA_interaction, leesedwards_distance<ndim> >(
            std::make_shared<sf_HS_WCA_interaction>(eps, sca),
            std::make_shared<leesedwards_distance<ndim> >(boxvec, shear),
            boxvec,
            2 * (1 + sca) * *std::max_element(radii.begin(), radii.end()), // rcut
            ncellx_scale, radii, sca, balance_omp)
    {}
    size_t get_nr_unique_pairs() const { return CellListPotential< sf_HS_WCA_interaction, leesedwards_distance<ndim> >::m_celliter->get_nr_unique_pairs(); }
};

/**
 * Pairwise WCA potential with interaction lists
 */
class HS_WCANeighborList : public SimplePairwiseNeighborList< sf_HS_WCA_interaction > {
public:
    HS_WCANeighborList(Array<size_t> & ilist, double eps, double sca, Array<double> radii)
        :  SimplePairwiseNeighborList< sf_HS_WCA_interaction > (
                std::make_shared<sf_HS_WCA_interaction>(eps, sca), ilist, radii)
    {
    }
};


} //namespace pele
#endif //#ifndef _PELE_HS_WCA_H
