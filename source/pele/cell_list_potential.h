#ifndef _PELE_CELL_LIST_POTENTIAL_H
#define _PELE_CELL_LIST_POTENTIAL_H

#include <iostream>
#include <memory>
#include <vector>
#include <utility>
#include <stdexcept>

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
    virtual double get_energy_gradient(Array<double> xa, Array<double> grad)
    {
        refresh_cell_list(xa);
        const double* x = xa.data();
        double result = 0;
        grad.assign(double(0));
        if (xa.size() != grad.size()) {
            throw std::runtime_error("CellListPotential::get_energy_gradient: illegal input");
        }
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
            double gij;
            result += m_interaction->energy_gradient(r2, &gij, i, j);
            for (size_t k = 0; k < m_ndim; ++k) {
                grad[xi_off + k] -= gij * dr[k];
            }
            for (size_t k = 0; k < m_ndim; ++k) {
                grad[xj_off + k] += gij * dr[k];
            }
            delete[] dr;
        }
        return result;
    }
    virtual double get_energy_gradient_hessian(Array<double> xa,
            Array<double> grad, Array<double> hess)
    {
        if (xa.size() != grad.size()) {
            throw std::runtime_error("CellListPotential::get_energy_gradient_hessian: illegal input grad");
        }
        if (hess.size() != xa.size() * xa.size()) {
            throw std::runtime_error("CellListPotential::get_energy_gradient_hessian: illegal input hess");
        }
        refresh_cell_list(xa);
        const double* x = xa.data();
        double result = 0;
        grad.assign(double(0));
        hess.assign(double(0));
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
            double gij;
            double hij;
            result += m_interaction->energy_gradient_hessian(r2, &gij, &hij, i, j);
            for (size_t k = 0; k < m_ndim; ++k) {
                grad[xi_off + k] -= gij * dr[k];
            }
            for (size_t k = 0; k < m_ndim; ++k) {
                grad[xj_off + k] += gij * dr[k];
            }
            //this part is copied from simple_pairwise_potential.h
            //(even more so than the rest)
            const size_t N = xa.size();
            const size_t i1 = xi_off;
            const size_t j1 = xj_off;
            for (size_t k = 0; k < m_ndim; ++k) {
                //diagonal block - diagonal terms
                const double Hii_diag = (hij + gij) * dr[k] * dr[k] / r2 - gij;
                hess[N * (i1 + k) + i1 + k] += Hii_diag;
                hess[N * (j1 + k) + j1 + k] += Hii_diag;
                //off diagonal block - diagonal terms
                const double Hij_diag = -Hii_diag;
                hess[N * (i1 + k) + j1 + k] = Hij_diag;
                hess[N * (j1 + k) + i1 + k] = Hij_diag;
                for (size_t l = k + 1; l < m_ndim; ++l) {
                    //diagonal block - off diagonal terms
                    const double Hii_off = (hij + gij) * dr[k] * dr[l] / r2;
                    hess[N * (i1 + k) + i1 + l] += Hii_off;
                    hess[N * (i1 + l) + i1 + k] += Hii_off;
                    hess[N * (j1 + k) + j1 + l] += Hii_off;
                    hess[N * (j1 + l) + j1 + k] += Hii_off;
                    //off diagonal block - off diagonal terms
                    const double Hij_off = -Hii_off;
                    hess[N * (i1 + k) + j1 + l] = Hij_off;
                    hess[N * (i1 + l) + j1 + k] = Hij_off;
                    hess[N * (j1 + k) + i1 + l] = Hij_off;
                    hess[N * (j1 + l) + i1 + k] = Hij_off;
                }
            }
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
