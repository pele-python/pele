#ifndef _PELE_FROZEN_ATOMS_H
#define _PELE_FROZEN_ATOMS_H

#include <vector>
#include <algorithm>
#include <set>
#include <assert.h>
#include <iostream>
#include <memory>

#include "array.h"
#include "base_potential.h"

namespace pele{
/**
 * Class for converting to and from a reduced representation and back.
 * The full set of coordinates includes the mobile and frozen degrees of
 * freedom.  The reduced set of coordinates includes only the mobile
 * degrees of freedom.
 */
class FrozenCoordsConverter
{
protected:
    std::vector<double> const _reference_coords;
    std::vector<size_t> _frozen_dof;
    std::vector<size_t> _mobile_dof;
public:
    FrozenCoordsConverter(Array<double> const reference_coords,
            Array<size_t> const frozen_dof) :
        _reference_coords(reference_coords.begin(), reference_coords.end())
    {
        //populate _frozen_dof after removing duplicates and sorting
        std::set<size_t> frozenset(frozen_dof.begin(), frozen_dof.end());
        _frozen_dof = std::vector<size_t>(frozenset.begin(), frozenset.end());
        std::sort(_frozen_dof.begin(), _frozen_dof.end());

        // do a sanity check
        if (_frozen_dof.size() > 0){
            if (_frozen_dof[_frozen_dof.size()-1] >= (size_t)ndof()) {
                throw std::runtime_error("index of frozen degree of freedom is out of bounds");
            }
        }

        //populate _mobile_dof
        _mobile_dof = std::vector<size_t>(ndof() - ndof_frozen());
        size_t imobile = 0;
        for (size_t i=0; i<_reference_coords.size(); ++i){
            // if degree of freedom i is not in frozen, add it to _mobile_dof
            if (frozenset.count(i) == 0){
                _mobile_dof[imobile] = i;
                ++imobile;
            }
        }
        assert(imobile == _mobile_dof.size());

    }

    size_t ndof() const { return _reference_coords.size(); }
    size_t ndof_frozen() const { return _frozen_dof.size(); }
    size_t ndof_mobile() const { return _mobile_dof.size(); }

    pele::Array<size_t> get_frozen_dof() { return pele::Array<size_t>(_frozen_dof); }
    pele::Array<size_t> get_mobile_dof() { return pele::Array<size_t>(_mobile_dof); }

    /**
     * Return the reduced representation of the system.  i.e. return a
     * set of coordinates with the frozen degrees of freedom removed.
     */
    Array<double> get_reduced_coords(Array<double> const &full_coords){
        if (full_coords.size() != ndof()) {
            std::invalid_argument("full_coords has the wrong size");
        }
        Array<double> reduced_coords(ndof_mobile());
        for (size_t i=0; i < ndof_mobile(); ++i){
            reduced_coords[i] = full_coords[_mobile_dof[i]];
        }
        return reduced_coords;
    }

    /**
     * Return the full representation of the system.  i.e. return a set
     * of coordinates with the frozen degrees of freedom added back in.
     */
    Array<double> get_full_coords(Array<double> const &reduced_coords)
    {
        if (reduced_coords.size() != ndof_mobile()) {
            std::invalid_argument("reduced_coords has the wrong size");
        }
        // make a new array full_coords as a copy of _reference_coords.
        //Array<double> const a(_reference_coords); //wrap _reference_coords in an Array (returns error due to _reference_coords being const)
        ////Array<double> full_coords(a.copy());
        Array<double> full_coords(ndof());
        std::copy(_reference_coords.begin(), _reference_coords.end(), full_coords.begin());
        // replace the mobile degrees of freedom with those in reduced_coords
        for (size_t i=0; i < _mobile_dof.size(); ++i){
            full_coords[_mobile_dof[i]] = reduced_coords[i];
        }
        return full_coords;
    }

    /**
     * Return the gradient of the full representation of the system.
     * The gradient of the frozen degrees of freedom will be set to zero.
     */
    Array<double> get_full_grad(Array<double> const &reduced_grad)
    {
        if (reduced_grad.size() != ndof_mobile()) {
            std::invalid_argument("reduced_grad has the wrong size");
        }
        Array<double> full_grad(ndof(), 0.);
        // replace the mobile degrees of freedom with those in reduced_grad
        for (size_t i=0; i < _mobile_dof.size(); ++i){
            full_grad[_mobile_dof[i]] = reduced_grad[i];
        }
        return full_grad;
    }

    /**
     * Return the Hessian of the reduced representation of the system.
     */
    Array<double> get_reduced_hessian(Array<double> const full_hess)
    {
        if (full_hess.size() != ndof()*ndof()) {
            std::invalid_argument("full_hessian has the wrong size");
        }
        Array<double> reduced_hess(ndof_mobile() * ndof_mobile());
//                size_t const N = ndof_mobile();
        for (size_t i=0; i < ndof_mobile(); ++i){
            for (size_t j=0; j < ndof_mobile(); ++j){
                size_t k_red = i * ndof_mobile() + j;
                size_t k_full = _mobile_dof[i] * ndof() + _mobile_dof[j];
                reduced_hess[k_red] = full_hess[k_full];
            }
        }
        return reduced_hess;
    }

};


class FrozenPotentialWrapper : public BasePotential
{
protected:
    std::shared_ptr<BasePotential> _underlying_potential;
public:
    FrozenCoordsConverter coords_converter;

    FrozenPotentialWrapper(std::shared_ptr<BasePotential> potential,
            Array<double> const reference_coords,
            Array<size_t> const frozen_dof) 
        : _underlying_potential(potential),
          coords_converter(reference_coords, frozen_dof)
          
    {}

    ~FrozenPotentialWrapper() {}

//    inline size_t ndof() const { return coords_converter.ndof(); }
//    inline size_t ndof_frozen() const { return coords_converter.ndof_frozen(); }
//    inline size_t ndof_mobile() const { return coords_converter.ndof_mobile(); }
    inline pele::Array<size_t> get_mobile_dof() { return coords_converter.get_mobile_dof(); }
    inline pele::Array<size_t> get_frozen_dof() { return coords_converter.get_frozen_dof(); }

    inline Array<double> get_reduced_coords(Array<double> const full_coords){
        return coords_converter.get_reduced_coords(full_coords);
    }
    Array<double> get_full_coords(Array<double> const reduced_coords){
        return coords_converter.get_full_coords(reduced_coords);
    }

    inline double get_energy(Array<double> reduced_coords) 
    {
        if (reduced_coords.size() != coords_converter.ndof_mobile()){
            throw std::runtime_error("reduced coords does not have the right size");
        }
        Array<double> full_coords(coords_converter.get_full_coords(reduced_coords));
        double const energy = _underlying_potential->get_energy(full_coords);
        //reduced_coords.assign(coords_converter.get_reduced_coords(full_coords));
        return energy;
    }

    inline double get_energy_gradient(Array<double> reduced_coords, Array<double> reduced_grad) 
    {
        if (reduced_coords.size() != coords_converter.ndof_mobile()){
            throw std::runtime_error("reduced coords does not have the right size");
        }
        if (reduced_grad.size() != coords_converter.ndof_mobile()) {
            throw std::invalid_argument("reduced_grad has the wrong size");
        }
        
        Array<double> full_coords(coords_converter.get_full_coords(reduced_coords));
        Array<double> gfull(coords_converter.ndof());
        double const energy = _underlying_potential->get_energy_gradient(full_coords, gfull);
        //reduced_coords.assign(coords_converter.get_reduced_coords(full_coords));
        reduced_grad.assign(coords_converter.get_reduced_coords(gfull));
        return energy;
    }

    inline double get_energy_gradient_hessian(Array<double> reduced_coords,
            Array<double> reduced_grad, Array<double> reduced_hess)
    {
        if (reduced_coords.size() != coords_converter.ndof_mobile()){
            throw std::runtime_error("reduced coords does not have the right size");
        }
        if (reduced_grad.size() != coords_converter.ndof_mobile()) {
            throw std::invalid_argument("reduced_grad has the wrong size");
        }
        if (reduced_hess.size() != coords_converter.ndof_mobile()*coords_converter.ndof_mobile()){
            throw std::invalid_argument("reduced_hess has the wrong size");
        }
        Array<double> full_coords(coords_converter.get_full_coords(reduced_coords));
        Array<double> gfull(coords_converter.ndof());
        Array<double> hfull(coords_converter.ndof()*coords_converter.ndof());
        const double energy = _underlying_potential->get_energy_gradient_hessian(full_coords, gfull, hfull);
        //reduced_coords.assign(coords_converter.get_reduced_coords(full_coords));
        reduced_grad.assign(coords_converter.get_reduced_coords(gfull));
        reduced_hess.assign(coords_converter.get_reduced_hessian(hfull));
        return energy;
    }
};
}

#endif
