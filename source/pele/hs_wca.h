#ifndef _PELE_HS_WCA_H
#define _PELE_HS_WCA_H

#include <algorithm>
#include <memory>

#include "simple_pairwise_potential.h"
#include "simple_pairwise_ilist.h"
#include "atomlist_potential.h"
#include "distance.h"
#include "frozen_atoms.h"
#include "cell_list_potential.h"

namespace pele {

/**
 * Fast pairwise interaction for Hard Sphere + Weeks-Chandler-Andersen (fHS_WCA) potential, refer to S. Martiniani CPGS pp 20
 * well depth _eps and scaling factor (shell thickness = sca * R, where R is the hard core radius),
 * sca determines the thickness of the shell
 */
struct HS_WCA_interaction {
    const double _eps;
    const double _sca;
    const double _infty;
    const double _prfac;
    const Array<double> _radii;

    HS_WCA_interaction(const double eps, const double sca, const Array<double> radii) 
        : _eps(eps),
          _sca(sca),
          _infty(std::pow(10.0, 50)),
          _prfac(std::pow((2 * _sca + _sca * _sca), 3) / std::sqrt(2)),
          _radii(radii.copy())
    {}

    /* calculate energy from distance squared, r0 is the hard core distance, r is the distance between the centres */
    double inline energy(const double r2, const size_t atomi, const size_t atomj) const 
    {
        const double r0 = _radii[atomi] + _radii[atomj]; //sum of the hard core radii
        const double r02 = r0 * r0;
        const double dr = r2 - r02; // note that dr is the difference of the squares
        const double ir2 = 1.0 / (dr * dr);
        const double ir6 = ir2 * ir2 * ir2;
        const double ir12 = ir6 * ir6;
        const double C3 = _prfac * r02 * r02 * r02;
        const double C6 = C3 * C3;
        const double C12 = C6 * C6;
        const double coff = r0 * (1.0 +_sca); //distance at which the soft cores are at contact
        if (r2 <= r02){
            return _infty;
        }
        if (r2 > coff * coff){
            return 0;
        }
        return 4. * _eps * (-C6 * ir6 + C12 * ir12) + _eps;
    }

    /* calculate energy and gradient from distance squared, gradient is in g/|rij|, r0 is the hard core distance, r is the distance between the centres */
    double inline energy_gradient(const double r2, double*const& gij, const size_t atomi, const size_t atomj) const 
    {
        const double r0 = _radii[atomi] + _radii[atomj]; //sum of the hard core radii
        const double r02 = r0 * r0;
        const double dr = r2 - r02; // note that dr is the difference of the squares
        const double ir2 = 1.0 / (dr * dr);
        const double ir6 = ir2 * ir2 * ir2;
        const double ir12 = ir6 * ir6;
        const double C3 = _prfac * r02 * r02 * r02;
        const double C6 = C3 * C3;
        const double C12 = C6 * C6;
        const double coff = r0 * (1.0 + _sca); //distance at which the soft cores are at contact
        if (r2 <= r02){
            *gij = _infty;
            return _infty;
        }
        if (r2 > coff * coff){
            *gij = 0.;
            return 0.;
        }
        *gij = _eps * (- 48. * C6 * ir6 + 96. * C12 * ir12) / dr; //this is -g/r, 1/dr because powers must be 7 and 13
        return 4. * _eps * (-C6 * ir6 + C12 * ir12) + _eps;
    }

    double inline energy_gradient_hessian(const double r2, double*const& gij, double*const& hij, const size_t atomi, const size_t atomj) const
    {
        const double r0 = _radii[atomi] + _radii[atomj]; //sum of the hard core radii
        const double r02 = r0 * r0;
        const double dr = r2 - r02; // note that dr is the difference of the squares
        const double ir2 = 1.0 / (dr * dr);
        const double ir6 = ir2 * ir2 * ir2;
        const double ir12 = ir6 * ir6;
        const double C3 = _prfac * r02 * r02 * r02;
        const double C6 = C3 * C3;
        const double C12 = C6 * C6;
        const double coff = r0 * (1.0 + _sca); //distance at which the soft cores are at contact
        if (r2 <= r02){
            *gij = _infty;
            *hij = _infty;
            return _infty;
        }
        if (r2 > coff * coff){
            *gij = 0.;
            *hij = 0.;
            return 0.;
        }
        *gij = _eps * (- 48. * C6 * ir6 + 96. * C12 * ir12) / dr; //this is -g/r, 1/dr because powers must be 7 and 13
        *hij = -*gij + _eps * ( -672. * C6 * ir6 + 2496. * C12 * ir12)  * r2 * ir2;
        return 4. * _eps * (-C6 * ir6 + C12 * ir12) + _eps;
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
class HS_WCA : public SimplePairwisePotential< HS_WCA_interaction, cartesian_distance<ndim> > {
public:
    HS_WCA(double eps, double sca, Array<double> radii)
        : SimplePairwisePotential< HS_WCA_interaction, cartesian_distance<ndim> >(
                std::make_shared<HS_WCA_interaction>(eps, sca, radii),
                std::make_shared<cartesian_distance<ndim> >()
            )
    {
        static_assert(ndim > 0, "illegal box dimension");
        if (eps < 0) {
            throw std::runtime_error("HS_WCA: illegal input: eps");
        }
        if (sca < 0) {
            throw std::runtime_error("HS_WCA: illegal inout: sca");
        }
        if (radii.size() == 0) {
            throw std::runtime_error("HS_WCA: illegal input: radii");
        }
    }
};

/**
 * Pairwise HS_WCA potential in a rectangular box
 */
template<size_t ndim>
class HS_WCAPeriodic : public SimplePairwisePotential< HS_WCA_interaction, periodic_distance<ndim> > {
public:
    HS_WCAPeriodic(double eps, double sca, Array<double> radii, Array<double> const boxvec)
        : SimplePairwisePotential< HS_WCA_interaction, periodic_distance<ndim> > (
                std::make_shared<HS_WCA_interaction>(eps, sca, radii),
                std::make_shared<periodic_distance<ndim> >(boxvec)
                )
    {
        static_assert(ndim > 0, "illegal box dimension");
        if (eps < 0) {
            throw std::runtime_error("HS_WCA: illegal input: eps");
        }
        if (sca < 0) {
            throw std::runtime_error("HS_WCA: illegal inout: sca");
        }
        if (radii.size() == 0) {
            throw std::runtime_error("HS_WCA: illegal input: radii");
        }
        if (boxvec.size() != ndim) {
            throw std::runtime_error("HS_WCA: illegal input: boxvec");
        }
    }
};

template<size_t ndim>
class HS_WCACellLists : public CellListPotential< HS_WCA_interaction, cartesian_distance<ndim> > {
public:
    HS_WCACellLists(double eps, double sca, Array<double> radii, Array<double> const boxvec,
            const double rcut, const double ncellx_scale = 1.0)
    : CellListPotential< HS_WCA_interaction, cartesian_distance<ndim> >(
            std::make_shared<HS_WCA_interaction>(eps, sca, radii),
            std::make_shared<cartesian_distance<ndim> >(),
            boxvec, rcut, ncellx_scale)
    {
        static_assert(ndim > 0, "illegal box dimension");
        if (eps < 0) {
            throw std::runtime_error("HS_WCA: illegal input: eps");
        }
        if (sca < 0) {
            throw std::runtime_error("HS_WCA: illegal inout: sca");
        }
        if (radii.size() == 0) {
            throw std::runtime_error("HS_WCA: illegal input: radii");
        }
        if (boxvec.size() != ndim) {
            throw std::runtime_error("HS_WCA: illegal input: boxvec");
        }
        if (rcut < 2 * (1 + sca) * *std::max_element(radii.data(), radii.data() + radii.size())) {
            throw std::runtime_error("HS_WCA: illegal input: rcut");
        }
    }
    size_t get_nr_unique_pairs() const { return CellListPotential< HS_WCA_interaction, cartesian_distance<ndim> >::m_celliter->get_nr_unique_pairs(); }
};

template<size_t ndim>
class HS_WCAPeriodicCellLists : public CellListPotential< HS_WCA_interaction, periodic_distance<ndim> > {
public:
    HS_WCAPeriodicCellLists(double eps, double sca, Array<double> radii, Array<double> const boxvec,
            const double rcut, const double ncellx_scale = 1.0)
    : CellListPotential< HS_WCA_interaction, periodic_distance<ndim> >(
            std::make_shared<HS_WCA_interaction>(eps, sca, radii),
            std::make_shared<periodic_distance<ndim> >(boxvec),
            boxvec, rcut, ncellx_scale)
    {
        static_assert(ndim > 0, "illegal box dimension");
        if (eps < 0) {
            throw std::runtime_error("HS_WCA: illegal input: eps");
        }
        if (sca < 0) {
            throw std::runtime_error("HS_WCA: illegal inout: sca");
        }
        if (radii.size() == 0) {
            throw std::runtime_error("HS_WCA: illegal input: radii");
        }
        if (boxvec.size() != ndim) {
            throw std::runtime_error("HS_WCA: illegal input: boxvec");
        }
        if (rcut < 2 * (1 + sca) * *std::max_element(radii.data(), radii.data() + radii.size())) {
            throw std::runtime_error("HS_WCA: illegal input: rcut");
        }
    }
    size_t get_nr_unique_pairs() const { return CellListPotential< HS_WCA_interaction, periodic_distance<ndim> >::m_celliter->get_nr_unique_pairs(); }
};

/**
 * Frozen particle HS_WCA potential
 */
template<size_t ndim>
class HS_WCAFrozen : public FrozenPotentialWrapper<HS_WCA<ndim> > {
public:
    HS_WCAFrozen(double eps, double sca, Array<double> radii, Array<double>& reference_coords, Array<size_t>& frozen_dof)
        : FrozenPotentialWrapper< HS_WCA<ndim> > ( std::make_shared<HS_WCA<ndim> >(eps, sca,
                    radii), reference_coords.copy(), frozen_dof.copy())
    {
        static_assert(ndim > 0, "illegal box dimension");
        if (eps < 0) {
            throw std::runtime_error("HS_WCA: illegal input: eps");
        }
        if (sca < 0) {
            throw std::runtime_error("HS_WCA: illegal inout: sca");
        }
        if (radii.size() == 0) {
            throw std::runtime_error("HS_WCA: illegal input: radii");
        }
        if (reference_coords.size() != ndim * radii.size()) {
            throw std::runtime_error("HS_WCA: illegal input: coords vs. radii");
        }
    }
};

/**
 * Frozen particle HS_WCAPeriodic potential
 */
template<size_t ndim>
class HS_WCAPeriodicFrozen : public FrozenPotentialWrapper<HS_WCAPeriodic<ndim> > {
public:
    HS_WCAPeriodicFrozen(double eps, double sca, Array<double> radii, 
            Array<double> const boxvec, Array<double>& reference_coords,
            Array<size_t>& frozen_dof)
        : FrozenPotentialWrapper< HS_WCAPeriodic<ndim> > (
                std::make_shared<HS_WCAPeriodic<ndim> >(eps, sca, radii, boxvec),
                reference_coords.copy(), frozen_dof.copy())
    {
        static_assert(ndim > 0, "illegal box dimension");
        if (eps < 0) {
            throw std::runtime_error("HS_WCA: illegal input: eps");
        }
        if (sca < 0) {
            throw std::runtime_error("HS_WCA: illegal inout: sca");
        }
        if (radii.size() == 0) {
            throw std::runtime_error("HS_WCA: illegal input: radii");
        }
        if (boxvec.size() != ndim) {
            throw std::runtime_error("HS_WCA: illegal input: boxvec");
        }
        if (reference_coords.size() != ndim * radii.size()) {
            throw std::runtime_error("HS_WCA: illegal input: coords vs. radii");
        }
    }
};

template<size_t ndim>
class HS_WCACellListsFrozen : public FrozenPotentialWrapper<HS_WCACellLists<ndim> > {
public:
    HS_WCACellListsFrozen(double eps, double sca, Array<double> radii,
            Array<double> const boxvec, Array<double>& reference_coords,
            Array<size_t>& frozen_dof, const double rcut, const double ncellx_scale = 1.0)
        : FrozenPotentialWrapper< HS_WCACellLists<ndim> > (
                std::make_shared<HS_WCACellLists<ndim> >(eps, sca, radii, boxvec, rcut, ncellx_scale),
                reference_coords.copy(), frozen_dof.copy())
    {
        static_assert(ndim > 0, "illegal box dimension");
        if (eps < 0) {
            throw std::runtime_error("HS_WCA: illegal input: eps");
        }
        if (sca < 0) {
            throw std::runtime_error("HS_WCA: illegal inout: sca");
        }
        if (radii.size() == 0) {
            throw std::runtime_error("HS_WCA: illegal input: radii");
        }
        if (boxvec.size() != ndim) {
            throw std::runtime_error("HS_WCA: illegal input: boxvec");
        }
        if (reference_coords.size() != ndim * radii.size()) {
            throw std::runtime_error("HS_WCA: illegal input: coords vs. radii");
        }
        if (rcut < 2 * (1 + sca) * *std::max_element(radii.data(), radii.data() + radii.size())) {
            throw std::runtime_error("HS_WCA: illegal input: rcut");
        }
    }
};

template<size_t ndim>
class HS_WCAPeriodicCellListsFrozen : public FrozenPotentialWrapper<HS_WCAPeriodicCellLists<ndim> > {
public:
    HS_WCAPeriodicCellListsFrozen(double eps, double sca, Array<double> radii,
            Array<double> const boxvec, Array<double>& reference_coords,
            Array<size_t>& frozen_dof, const double rcut, const double ncellx_scale = 1.0)
        : FrozenPotentialWrapper< HS_WCAPeriodicCellLists<ndim> > (
                std::make_shared<HS_WCAPeriodicCellLists<ndim> >(eps, sca, radii, boxvec, rcut, ncellx_scale),
                reference_coords.copy(), frozen_dof.copy())
    {
        static_assert(ndim > 0, "illegal box dimension");
        if (eps < 0) {
            throw std::runtime_error("HS_WCA: illegal input: eps");
        }
        if (sca < 0) {
            throw std::runtime_error("HS_WCA: illegal inout: sca");
        }
        if (radii.size() == 0) {
            throw std::runtime_error("HS_WCA: illegal input: radii");
        }
        if (boxvec.size() != ndim) {
            throw std::runtime_error("HS_WCA: illegal input: boxvec");
        }
        if (reference_coords.size() != ndim * radii.size()) {
            throw std::runtime_error("HS_WCA: illegal input: coords vs. radii");
        }
        if (rcut < 2 * (1 + sca) * *std::max_element(radii.data(), radii.data() + radii.size())) {
            throw std::runtime_error("HS_WCA: illegal input: rcut");
        }
    }
};

/**
 * Pairwise WCA potential with interaction lists
 */
class HS_WCANeighborList : public SimplePairwiseNeighborList< HS_WCA_interaction > {
public:
    HS_WCANeighborList(Array<size_t> & ilist, double eps, double sca, Array<double> radii)
        :  SimplePairwiseNeighborList< HS_WCA_interaction > (
                std::make_shared<HS_WCA_interaction>(eps, sca, radii), ilist)
    {
        if (eps < 0) {
            throw std::runtime_error("HS_WCA: illegal input: eps");
        }
        if (sca < 0) {
            throw std::runtime_error("HS_WCA: illegal inout: sca");
        }
        if (radii.size() == 0) {
            throw std::runtime_error("HS_WCA: illegal input: radii");
        }
    }
};

} //namespace pele
#endif //#ifndef _PELE_HS_WCA_H
