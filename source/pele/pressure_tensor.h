#ifndef _PELE_PRESSURE_TENSOR_H
#define _PELE_PRESSURE_TENSOR_H

#include "pele/base_potential.h"

namespace pele {

double pressure_tensor(std::shared_ptr<pele::BasePotential> pot,
                       pele::Array<double> x,
                       pele::Array<double> ptensor,
                       const double volume);
    
} // namespace pele

#endif // #ifndef _PELE_PRESSURE_TENSOR_H
