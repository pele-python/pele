#ifndef PYGMIN_POTENTIAL_H
#define PYGMIN_POTENTIAL_H

#include <assert.h>
#include <vector>
#include <stdexcept>
#include <iostream>
#include "array.h"

namespace pele {
    /***
     * basic potential interface for native potentials
     */
    class BasePotential {
    public:
        virtual ~BasePotential() {}

        /**
         * Return the energy of configuration x.  This is the only function which
         * must be overloaded
         */
        virtual double get_energy(Array<double> x)
        {
            throw std::runtime_error("BasePotential::get_energy must be overloaded");
        }

        /**
         * compute the energy and gradient, but don't intialize the gradient to zero
         */
        virtual double add_energy_gradient(Array<double> x, Array<double> grad)
        {
            throw std::runtime_error("BasePotential::add_energy_gradient must be overloaded");
        }

        /**
         * compute the energy and gradient.  If not overloaded it will compute the numerical gradient
         */
        virtual double get_energy_gradient(Array<double> x, Array<double> grad)
        {
            double energy = get_energy(x);
            numerical_gradient(x, grad);
            return energy;
        }

        /**
         * compute the numerical gradient
         */
        virtual void numerical_gradient(Array<double> x, Array<double> grad, double eps=1e-6)
        {
            if (x.size() != grad.size()) {
                throw std::invalid_argument("grad.size() be the same as x.size()");
            }

            Array<double> xnew(x.copy());
            for (size_t i=0; i<xnew.size(); ++i){
                xnew[i] -= eps;
                double eminus = get_energy(xnew);
                xnew[i] += 2. * eps;
                double eplus = get_energy(xnew);
                grad[i] = (eplus - eminus) / (2. * eps);
                xnew[i] = x[i];
            }
        }

        /**
         * compute the numerical gradient
         */
        virtual void numerical_hessian(Array<double> x, Array<double> hess, double eps=1e-6)
        {
            if (hess.size() != x.size()*x.size()) {
                throw std::invalid_argument("hess.size() be the same as x.size()");
            }
            size_t const N = x.size();

            Array<double> gplus(x.size());
            Array<double> gminus(x.size());

            for (size_t i=0; i<x.size(); ++i){
                double xbackup = x[i];
                x[i] -= eps;
                get_energy_gradient(x, gminus);
                x[i] += 2. * eps;
                get_energy_gradient(x, gplus);
                x[i] = xbackup;

                for (size_t j=0; j<x.size(); ++j){
                    // hess[i,:] = (g1 - g2) / (2. * eps)
                    hess[N*i + j] = (gplus[j] - gminus[j]) / (2.*eps);
//                    hess[N*i + j] = (gplus[j] - gminus[j]) / (2.*eps);
                }
            }
        }

    };
}

#endif
