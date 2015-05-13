#include "pele/pressure_tensor.h"

namespace pele {

/*
 * computes the static pressure tensor, ignoring the momenta of the atoms
 * the momentum component can be added
 * */
double pressure_tensor(std::shared_ptr<pele::BasePotential> pot,
                       pele::Array<double> x,
                       pele::Array<double> ptensor,
                       const double volume)
{
    const size_t ndim = pot->get_ndim();
    const size_t natoms = x.size() / ndim;
    if (ndim * natoms != x.size()) {
        throw std::runtime_error("x is not divisible by the number of dimensions");
    }
    if (ptensor.size() != ndim * ndim) {
        throw std::runtime_error("ptensor must have size ndim * ndim");
    }

    double gij;
    double dr[ndim];
    ptensor.assign(0.);

    for (size_t atomi = 0; atomi < natoms; ++atomi) {
        size_t const i1 = ndim * atomi;
        for (size_t atomj = 0; atomj < atomi; ++atomj) {
            size_t const j1 = ndim * atomj;

            pot->get_rij(dr, &x[i1], &x[j1]);

            double r2 = 0;
            for (size_t k = 0; k < ndim; ++k) {
                r2 += dr[k] * dr[k];
            }
            double e = pot->get_interaction_energy_gradient(r2, &gij, atomi, atomj);

            for (size_t k = 0; k < ndim; ++k) {
                for (size_t l = k; l < ndim; ++l) {
                    ptensor[k * ndim + l] -= dr[l] * gij * dr[k];
                    ptensor[l * ndim + k] -= dr[k] * gij * dr[l]; //pressure tensor is symmetric
                }
            }
        }
    }
    ptensor /= volume;

    //pressure is the average of the trace of the pressure tensor
    double traceP = 0.;
    for (size_t i = 0; i < ndim; ++i) {
        traceP += ptensor[i + ndim + i];
    }
    return traceP / ndim;
}
    
} // namespace pele
