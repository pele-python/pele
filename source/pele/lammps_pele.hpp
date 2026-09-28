#ifndef LAMMPS_PELE_H
#define LAMMPS_PELE_H

#include "pele/base_potential.h"
#include <cstdint>

namespace pele {

class LAMMPSPotential : public BasePotential {
protected:
    uintptr_t handle;

public:
    LAMMPSPotential(uintptr_t handle);
    double get_energy(Array<double> xs);
    double get_energy_gradient(Array<double> xs, Array<double> gs);
};

}

#endif