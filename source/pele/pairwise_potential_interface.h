#ifndef PYGMIN_PAIRWISE_POTENTIAL_INTERFACE_H
#define PYGMIN_PAIRWISE_POTENTIAL_INTERFACE_H

#include "base_potential.h"
#include "array.h"

namespace pele {

class PairwisePotentialInterface : public BasePotential {
public:
    virtual ~PairwisePotentialInterface() {}
    /**
     * Return the number of dimensions (box dimensions).
     * Ideally this should be overloaded.
     */
    virtual inline size_t get_ndim() const
    {
        throw std::runtime_error("PairwisePotentialInterface::get_ndim must be overloaded");
    }

    /**
     * Return the distance as measured by the distance policy.
     */
    virtual inline void get_rij(double * const r_ij, double const * const r1, double const * const r2) const
    {
        throw std::runtime_error("PairwisePotentialInterface::get_rij must be overloaded");
    }

    /**
     * Return energy_gradient of interaction.
     */
    virtual inline double get_interaction_energy_gradient(double r2, double *gij, size_t atom_i, size_t atom_j) const
    {
        throw std::runtime_error("PairwisePotentialInterface::get_interaction_energy_gradient must be overloaded");
    }

    /**
     * Return gradient and Hessian of interaction.
     */
    virtual inline double get_interaction_energy_gradient_hessian(double r2, double *gij, double *hij, size_t atom_i, size_t atom_j) const
    {
        throw std::runtime_error("PairwisePotentialInterface::get_interaction_energy_gradient_hessian must be overloaded");
    }

    /**
     * Return lists of neighbours considering only certain atoms.
     */
    virtual void get_neighbours_picky(pele::Array<double> & coords,
                                      pele::Array<std::vector<size_t>> & neighbour_indss,
                                      pele::Array<std::vector<std::vector<double>>> & neighbour_distss,
                                      pele::Array<short> const & include_atoms,
                                      const double cutoff_factor = 1.0)
    {
        throw std::runtime_error("PairwisePotentialInterface::get_neighbours_picky must be overloaded");
    }

    /**
     * Return lists of neighbours.
     */
    virtual void get_neighbours(pele::Array<double> & coords,
                                pele::Array<std::vector<size_t>> & neighbour_indss,
                                pele::Array<std::vector<std::vector<double>>> & neighbour_distss,
                                const double cutoff_factor = 1.0)
    {
        throw std::runtime_error("PairwisePotentialInterface::get_neighbours must be overloaded");
    }

    /**
     * Return linearized vector of overlapping atoms.
     *
     * Every two entries in the resulting vector correspond to a pair of overlapping atoms.
     */
    virtual std::vector<size_t> get_overlaps(Array<double> & coords)
    {
        throw std::runtime_error("PairwisePotentialInterface::get_overlaps must be overloaded");
    }
};

}

#endif // #ifndef PYGMIN_PAIRWISE_POTENTIAL_INTERFACE_H
