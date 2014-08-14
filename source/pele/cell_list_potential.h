#ifndef _PELE_CELL_LIST_POTENTIAL_H
#define _PELE_CELL_LIST_POTENTIAL_H

#include <iostream>
#include <memory>
#include <vector>
#include <utility>

#include "base_potential.h"
#include "array.h"
#include "distance.h"
#include "neighbor_iterator.h"

namespace pele{
/**
 * Potential to iterate over the list of atom pairs generated with the
 * cell list implementation in neighbor_iterator.h.
 * This should also do the cell list construction and refresh, such that
 * the iterface is the same for the user as with SimplePairwise.
 */
template <typename pairwise_interaction, typename distance_policy>
class CellListPotential : public BasePotential {
protected:
    const static size_t m_ndim = distance_policy::_ndim;
    std::shared_ptr<pairwise_interaction> m_interaction;
    std::shared_ptr<distance_policy> m_dist;
    std::shared_ptr<CellIter<distance_policy> > m_celliter;
public:
    virtual ~CellListPotential() {}
    CellListPotential(std::shared_ptr<pairwise_interaction> interaction,
            std::shared_ptr<distance_policy> dist,
            std::shared_ptr<CellIter<distance_policy> > celliter)
        : m_interaction(interaction),
          m_dist(dist),
          m_celliter(celliter)
    {}
    virtual double get_energy(Array<double> xa)
    {
        refresh_cell_list(xa);
        const double* x = xa.data();
        double result = 0;
        for (auto ijpair = m_celliter->begin(); ijpair != m_celliter->end(); ++ijpair) {
            const size_t i = ijpair->first;
            const size_t j = ijpair->second;
            const size_t xi_off = m_ndim * i;
            const size_t xj_off = m_ndim * j;
            double* dr = new double[m_ndim];
            m_dist->get_rij(dr, x + xi_off, x + xj_off);
            double r2 = 0;
            for (size_t k = 0; k < m_ndim; ++k) {
                r2 += dr[k] * dr[k];
            }
            result += m_interaction->energy(r2, i, j);
            delete[] dr;
        }
        return result;
    }
protected:
    void refresh_cell_list(Array<double> x)
    {
        m_celliter->reset(x);
    }
}; //class CellListPotential

} //namespace pele

#endif //#ifndef _PELE_CELL_LIST_POTENTIAL_H
