#ifndef PYGMIN_SIMPLE_PAIRWISE_POTENTIAL_H
#define PYGMIN_SIMPLE_PAIRWISE_POTENTIAL_H

#include "pairwise_potential_interface.h"
#include "array.h"
#include "distance.h"
#include <memory>

namespace pele {

/**
 * Define a base class for potentials with simple pairwise interactions that
 * depend only on magnitude of the atom separation
 *
 * This class loops though atom pairs, computes the distances and get's the
 * value of the energy and gradient from the class pairwise_interaction.
 * pairwise_interaction is a passed parameter and defines the actual
 * potential function.
 */
template<typename pairwise_interaction,
    typename distance_policy = cartesian_distance<3> >
class SimplePairwisePotential : public PairwisePotentialInterface
{
protected:
    static const size_t m_ndim = distance_policy::_ndim;
    std::shared_ptr<pairwise_interaction> _interaction;
    std::shared_ptr<distance_policy> _dist;
    const double m_radii_sca;

    SimplePairwisePotential( std::shared_ptr<pairwise_interaction> interaction,
            std::shared_ptr<distance_policy> dist=NULL,
            const double radii_sca=0.0)
        : _interaction(interaction),
          _dist(dist),
          m_radii_sca(radii_sca)
    {
        if(_dist == NULL) _dist = std::make_shared<distance_policy>();
    }

public:
    virtual ~SimplePairwisePotential()
    {}
    virtual inline size_t get_ndim() const { return m_ndim; }

    virtual double get_energy(Array<double> & x);
    virtual double get_energy_gradient(Array<double> & x, Array<double> & grad)
    {
        grad.assign(0);
        return add_energy_gradient(x, grad);
    }
    virtual double get_energy_gradient_hessian(Array<double> & x, Array<double> & grad, Array<double> & hess)
    {
        grad.assign(0);
        hess.assign(0);
        return add_energy_gradient_hessian(x, grad, hess);
    }
    virtual double add_energy_gradient(Array<double> & x, Array<double> & grad);
    virtual double add_energy_gradient_hessian(Array<double> & x, Array<double> & grad, Array<double> & hess);
    virtual void get_neighbours(Array<double> & coords,
                                pele::Array<std::vector<size_t>> & neighbour_indss,
                                pele::Array<std::vector<std::vector<double>>> & neighbour_distss);
    virtual inline void get_rij(double * const r_ij, double const * const r1, double const * const r2) const
    {
        return _dist->get_rij(r_ij, r1, r2);
    }
    virtual inline double get_interaction_energy_gradient(double r2, double *gij, size_t atom_i, size_t atom_j) const
    {
        return _interaction->energy_gradient(r2, gij, atom_i, atom_j);
    }
    virtual inline double get_interaction_energy_gradient_hessian(double r2, double *gij, double *hij, size_t atom_i, size_t atom_j) const
    {
        return _interaction->energy_gradient_hessian(r2, gij, hij, atom_i, atom_j);
    }
};

template<typename pairwise_interaction, typename distance_policy>
inline double
SimplePairwisePotential<pairwise_interaction,distance_policy>::add_energy_gradient(
        Array<double> & x, Array<double> & grad)
{
    const size_t natoms = x.size() / m_ndim;
    if (m_ndim * natoms != x.size()) {
        throw std::runtime_error("x is not divisible by the number of dimensions");
    }
    if (grad.size() != x.size()) {
        throw std::runtime_error("grad must have the same size as x");
    }

    double e = 0.;
    double gij;
    double dr[m_ndim];

//    grad.assign(0.);

    for (size_t atomi=0; atomi<natoms; ++atomi) {
        size_t const i1 = m_ndim * atomi;
        for (size_t atomj=0; atomj<atomi; ++atomj) {
            size_t const j1 = m_ndim * atomj;

            _dist->get_rij(dr, &x[i1], &x[j1]);

            double r2 = 0;
            for (size_t k=0; k<m_ndim; ++k) {
                r2 += dr[k]*dr[k];
            }
            e += _interaction->energy_gradient(r2, &gij, atomi, atomj);
            for (size_t k=0; k<m_ndim; ++k) {
                grad[i1+k] -= gij * dr[k];
            }
            for (size_t k=0; k<m_ndim; ++k) {
                grad[j1+k] += gij * dr[k];
            }
        }
    }
    return e;
}

template<typename pairwise_interaction, typename distance_policy>
inline double SimplePairwisePotential<pairwise_interaction, distance_policy>::add_energy_gradient_hessian(Array<double> & x,
        Array<double> & grad, Array<double> & hess)
{
    double hij, gij;
    double dr[m_ndim];
    const size_t N = x.size();
    const size_t natoms = x.size()/m_ndim;
    if (m_ndim * natoms != x.size()) {
        throw std::runtime_error("x is not divisible by the number of dimensions");
    }
    if (x.size() != grad.size()) {
        throw std::invalid_argument("the gradient has the wrong size");
    }
    if (hess.size() != x.size() * x.size()) {
        throw std::invalid_argument("the Hessian has the wrong size");
    }
//    hess.assign(0.);
//    grad.assign(0.);

    double e = 0.;
    for (size_t atomi=0; atomi<natoms; ++atomi) {
        int i1 = m_ndim*atomi;
        for (size_t atomj=0;atomj<atomi;++atomj){
            int j1 = m_ndim*atomj;
            _dist->get_rij(dr, &x[i1], &x[j1]);
            double r2 = 0;
            for (size_t k=0;k<m_ndim;++k){r2 += dr[k]*dr[k];}

            e += _interaction->energy_gradient_hessian(r2, &gij, &hij, atomi, atomj);

            for (size_t k=0; k<m_ndim; ++k)
                grad[i1+k] -= gij * dr[k];
            for (size_t k=0; k<m_ndim; ++k)
                grad[j1+k] += gij * dr[k];

            for (size_t k=0; k<m_ndim; ++k){
                //diagonal block - diagonal terms
                double Hii_diag = (hij+gij)*dr[k]*dr[k]/r2 - gij;
                hess[N*(i1+k)+i1+k] += Hii_diag;
                hess[N*(j1+k)+j1+k] += Hii_diag;
                //off diagonal block - diagonal terms
                double Hij_diag = -Hii_diag;
                hess[N*(i1+k)+j1+k] = Hij_diag;
                hess[N*(j1+k)+i1+k] = Hij_diag;
                for (size_t l = k+1; l<m_ndim; ++l){
                    //diagonal block - off diagonal terms
                    double Hii_off = (hij+gij)*dr[k]*dr[l]/r2;
                    hess[N*(i1+k)+i1+l] += Hii_off;
                    hess[N*(i1+l)+i1+k] += Hii_off;
                    hess[N*(j1+k)+j1+l] += Hii_off;
                    hess[N*(j1+l)+j1+k] += Hii_off;
                    //off diagonal block - off diagonal terms
                    double Hij_off = -Hii_off;
                    hess[N*(i1+k)+j1+l] = Hij_off;
                    hess[N*(i1+l)+j1+k] = Hij_off;
                    hess[N*(j1+k)+i1+l] = Hij_off;
                    hess[N*(j1+l)+i1+k] = Hij_off;
                }
            }
        }
    }
    return e;
}

template<typename pairwise_interaction, typename distance_policy>
inline double SimplePairwisePotential<pairwise_interaction, distance_policy>::get_energy(Array<double> & x)
{
    size_t const natoms = x.size()/m_ndim;
    if (m_ndim * natoms != x.size()) {
        throw std::runtime_error("x is not divisible by the number of dimensions");
    }
    double e=0.;
    double dr[m_ndim];

    for (size_t atomi=0; atomi<natoms; ++atomi) {
        size_t i1 = m_ndim*atomi;
        for (size_t atomj=0; atomj<atomi; ++atomj) {
            size_t j1 = m_ndim*atomj;
            _dist->get_rij(dr, &x[i1], &x[j1]);
            double r2 = 0;
            for (size_t k=0;k<m_ndim;++k) {
                r2 += dr[k]*dr[k];
            }
            e += _interaction->energy(r2, atomi, atomj);
        }
    }
    return e;
}

template<typename pairwise_interaction, typename distance_policy>
void SimplePairwisePotential<pairwise_interaction, distance_policy>::get_neighbours(
    Array<double> & coords,
    pele::Array<std::vector<size_t>> & neighbour_indss,
    pele::Array<std::vector<std::vector<double>>> & neighbour_distss)
{
    size_t natoms = coords.size()/m_ndim;
    if (m_ndim * natoms != coords.size()) {
        throw std::runtime_error("coords is not divisible by the number of dimensions");
    }
    if (_interaction->m_radii.size() == 0) {
        throw std::runtime_error("Can't calculate neighbours, because the "
                                 "used interaction doesn't use radii. ");
    }
    std::vector<double> dr(m_ndim);
    std::vector<double> neg_dr(m_ndim);
    neighbour_indss = pele::Array<std::vector<size_t>>(natoms);
    neighbour_distss = pele::Array<std::vector<std::vector<double>>>(natoms);

    for (size_t atomi=0; atomi<natoms; ++atomi) {
        size_t i1 = m_ndim*atomi;
        for (size_t atomj=0; atomj<atomi; ++atomj) {
            size_t j1 = m_ndim*atomj;
            _dist->get_rij(dr.data(), &coords[i1], &coords[j1]);
            double r2 = 0;
            for (size_t k=0;k<m_ndim;++k) {
                r2 += dr[k]*dr[k];
                neg_dr[k] = -dr[k];
            }
            const double r_H = _interaction->m_radii[atomi] + _interaction->m_radii[atomj];
            const double r_S = (1 + m_radii_sca) * r_H;
            const double r_S2 = r_S * r_S;
            if(r2 <= r_S2) {
                neighbour_indss[atomi].push_back(atomj);
                neighbour_indss[atomj].push_back(atomi);
                neighbour_distss[atomi].push_back(dr);
                neighbour_distss[atomj].push_back(neg_dr);
            }
        }
    }
}

} // namespace pele

#endif // #ifndef PYGMIN_SIMPLE_PAIRWISE_POTENTIAL_H
