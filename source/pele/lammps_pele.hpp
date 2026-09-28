#ifndef LAMMPS_PELE_H
#define LAMMPS_PELE_H

#include "pele/base_potential.hpp"
#include <cstdint>

namespace pele {

class LAMMPSPotential : public BasePotential {
protected:
    uintptr_t handle;

public:
    LAMMPSPotential(uintptr_t handle);
    using BasePotential::get_energy_gradient;
    double get_energy(Array<double> const &xs) override;
    double get_energy_gradient(Array<double> const &xs, Array<double> &gs) override;
};

}

#endif
