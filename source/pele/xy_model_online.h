#ifndef _PELE_XY_MODEL_ONLINE_H
#define _PELE_XY_MODEL_ONLINE_H

#include "pele/base_potential_online.h"

namespace pele {

class XYModelOnline : public BasePotentialOnline {
public:
    XYModelOnline(const size_t nr_spins)
        : BasePotentialOnline(nr_spins)
    {}
};
    
} // namespace pele

#endif // #ifndef _PELE_XY_MODEL_ONLINE_H

