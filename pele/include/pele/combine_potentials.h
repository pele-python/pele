#ifndef _PELE_COMBINE_POTENTIALS_H_
#define _PELE_COMBINE_POTENTIALS_H_
#include <list>
#include <memory>
#include "base_potential.h"

namespace pele {

/**
 * Potential wrapper which wraps multiple potentials to
 * act as one potential.  This can be used to implement
 * multiple types of interactions, e.g. a system with two
 * types atoms
 */
class CombinedPotential : public BasePotential{
protected:
    std::list<std::shared_ptr<BasePotential> > _potentials;

public:
    CombinedPotential() {}

    /**
     * destructor: destroy all the potentials in the list
     */
    virtual ~CombinedPotential() {}

    /**
     * add a potential to the list
     */
    virtual void add_potential(std::shared_ptr<BasePotential> potential)
    {
        _potentials.push_back(potential);
    }

    virtual double get_energy(Array<double> x)
    {
        double energy = 0.;
        for (auto & pot_ptr : _potentials){
            energy += pot_ptr->get_energy(x);
        }
        return energy;
    }

    virtual double get_energy_gradient(Array<double> x, Array<double> grad)
    {
        if (x.size() != grad.size()) {
            throw std::invalid_argument("the gradient has the wrong size");
        }

        double energy = 0.;
        grad.assign(0.);

        for (auto & pot_ptr : _potentials){
            energy += pot_ptr->add_energy_gradient(x, grad);
        }
        return energy;
    }

    virtual double get_energy_gradient_hessian(Array<double> x, Array<double> grad,
            Array<double> hess)
    {
        if (x.size() != grad.size()) {
            throw std::invalid_argument("the gradient has the wrong size");
        }
        if (hess.size() != x.size() * x.size()) {
            throw std::invalid_argument("the Hessian has the wrong size");
        }

        double energy = 0.;
        grad.assign(0.);
        hess.assign(0.);

        for (auto & pot_ptr : _potentials){
            energy += pot_ptr->add_energy_gradient_hessian(x, grad, hess);
        }
        return energy;
    }

};
}

#endif
