#ifndef _PELE_LOWEST_EIG_POTENTIAL_H
#define _PELE_LOWEST_EIG_POTENTIAL_H

#include "base_potential.h"
#include "array.h"
#include <algorithm>
#include <vector>
#include <cmath>
#include <memory>

namespace pele {

inline void zero_modes_translational(std::vector<pele::Array<double> > & zev,
        size_t natoms, size_t bdim)
{
    double v = 1 / sqrt(natoms);
    size_t N = natoms * bdim;
    for(size_t i=0; i<bdim; ++i) {
        pele::Array<double> evec(N, 0); //initialize array of zeros
        for(size_t j=i; j<N; j+=bdim) {
            evec[j] = v;
        }
        zev.push_back(evec);
    }
}

class Orthogonalize{
public:
    virtual ~Orthogonalize(){};
    virtual void orthogonalize(Array<double>& coords, Array<double>& vector) =0;
};

class OrthogonalizeTranslational : public Orthogonalize{
protected:
    std::vector<pele::Array<double>> _tr_evec;
    size_t _natoms, _bdim, _ndim;
    double _d, _tol;
public:

    //OrthogonalizeTranslational(size_t natoms, size_t bdim, double tol=1e-6);
    OrthogonalizeTranslational(size_t natoms, size_t bdim, double tol=1e-6)
        : _natoms(natoms), _bdim(bdim), _ndim(bdim*natoms), _tol(tol)
    {
        /*initialize translational eigenvectors to canonical orthonormal basis*/
        zero_modes_translational(_tr_evec, _natoms, _bdim);
    }

    virtual ~OrthogonalizeTranslational() {}

    //virtual inline void orthogonalize(Array<double>& coords, Array<double>& vector);
    virtual inline void orthogonalize(Array<double>& coords, Array<double>& vector)
    {
        bool success = true;
        pele::Array<double> dot_prod(_bdim);
        vector /= norm(vector);
        //generally in this loop success will be set to false
        for (size_t i=0; i<_bdim;++i) {
          dot_prod[i] = dot(_tr_evec[i],vector);
          if(std::abs(dot_prod[i]) > _tol){success = false;};
        }

        while (success == false) {
            success = true;
            for (size_t i=0; i<_bdim;++i) {
                for(size_t j=0;j<_ndim;++j) {
                    vector[j] -= dot_prod[i]*_tr_evec[i][j];
                }
            }
            vector /= norm(vector);
            for (size_t i=0; i<_bdim;++i) {
                dot_prod[i] = dot(_tr_evec[i],vector);
                if (std::abs(dot_prod[i]) > _tol) {
                    success = false;
                };
            }
        }
    }

}; //class OrthogonalizeTranslational

/*
 * Lowest Eigenvalue Potential:
 * gd = g(x+d)
 * m_v: normal unit vector in arbitrary direction
 * m_d: finite difference step
 * */

class LowestEigPotential : public BasePotential {
protected:
    std::shared_ptr<pele::BasePotential> m_potential;
    pele::Array<double> m_coords, m_coordsd_plus, m_coordsd_minus, m_g, m_gdplus, m_gdminus, m_gdiff;
    size_t m_bdim, m_natoms;
    double m_d;
    bool m_first_order;
    OrthogonalizeTranslational m_orthog;
public:

    //LowestEigPotential(std::shared_ptr<pele::BasePotential> potential, pele::Array<double>
    //        coords, size_t bdim, double d=1e-6);
    /*constructor*/
    LowestEigPotential(std::shared_ptr<pele::BasePotential> potential, pele::Array<double>coords,
            size_t bdim, double d=1e-6, bool first_order=true)
        : m_potential(potential),
          m_coords(coords),
          m_coordsd_plus(coords.size()),
          m_coordsd_minus(coords.size()),
          m_g(m_coords.size()),
          m_gdplus(m_coords.size()),
          m_gdminus(m_coords.size()),
          m_gdiff(m_coords.size()),
          m_bdim(bdim),
          m_natoms(m_coords.size()/m_bdim),
          m_d(d),
          m_first_order(first_order),
          m_orthog(m_natoms,m_bdim)
    {
        m_potential->get_energy_gradient(m_coords,m_g);
    }


    virtual ~LowestEigPotential(){}


    double m_curvature_first_order(pele::Array<double> x){
        for (size_t i=0;i<x.size();++i) {
            m_coordsd_plus[i] = m_coords[i] + m_d*x[i];
        }
        m_potential->get_energy_gradient(m_coordsd_plus, m_gdplus);
        for (size_t i=0;i<m_gdiff.size();++i) {
            m_gdiff[i] = m_gdplus[i] - m_g[i];
        }
        return dot(m_gdiff,x)/m_d;
    }

    double m_curvature_second_order(pele::Array<double> x){
        for (size_t i=0;i<x.size();++i) {
            double dx = m_d*x[i];
            m_coordsd_plus[i] = m_coords[i] + dx;
            m_coordsd_minus[i] = m_coords[i] - dx;
        }

        m_potential->get_energy_gradient(m_coordsd_plus,m_gdplus);
        m_potential->get_energy_gradient(m_coordsd_minus,m_gdminus);

        for (size_t i=0;i<x.size();++i) {
            m_gdiff[i] = m_gdplus[i] - m_gdminus[i];
        }
        return dot(m_gdiff,x)/(2. * m_d);
    }

    double m_curvature(pele::Array<double> x){
        if (m_first_order){
            return m_curvature_first_order(x);
        }
        else{
            return m_curvature_second_order(x);
        }
    }

    //virtual double inline get_energy(pele::Array<double> x);
    /* calculate energy from distance squared, r0 is the hard core distance, r is the distance between the centres */
    virtual double inline get_energy(pele::Array<double> x)
    {
        m_orthog.orthogonalize(m_coords, x); //takes care of orthogonalizing and normalizing x
        return this->m_curvature(x);
    }

    //virtual double inline get_energy_gradient(pele::Array<double> x,
    //        pele::Array<double> grad);
    /* calculate energy and gradient from distance squared, gradient is in g/|rij|, r0 is the hard core distance, r is the distance between the centres */
    virtual double inline get_energy_gradient(pele::Array<double> x, pele::Array<double> grad)
    {
        m_orthog.orthogonalize(m_coords, x);  //takes care of orthogonalizing and normalizing x

        double mu = this->m_curvature(x);

        if (m_first_order){
            for (size_t i=0;i<x.size();++i) {
                grad[i] = 2.*m_gdiff[i]/m_d - 2.*mu*x[i];
            }
        }
        else{
            for (size_t i=0;i<x.size();++i) {
                grad[i] = m_gdiff[i]/m_d - 2.*mu*x[i];
            }
        }

        return mu;
    }

    //void reset_coords(pele::Array<double> new_coords);
    void reset_coords(pele::Array<double> new_coords)
    {
        m_coords.assign(new_coords);
        m_potential->get_energy_gradient(m_coords, m_g);
    }


};

}//namespace pele

#endif//#ifndef _PELE_LOWEST_EIG_POTENTIAL_H
