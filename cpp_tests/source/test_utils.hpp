#ifndef __PELE_TEST_UTILS_HPP__
#define __PELE_TEST_UTILS_HPP__

#include "pele/array.h"
#include "pele/base_potential.h"

#include <iostream>
#include <stdexcept>
#include <gtest/gtest.h>
#include <cmath>
#include <memory>

using pele::Array;

class PotentialTest :  public ::testing::Test
{
public:
    std::shared_ptr<pele::BasePotential> pot; // this must be set manually
    pele::Array<double> x; //this must be set manually
    double etrue; // this must be set manually
    pele::Array<double> g, gnum, h, hnum;

    void test_energy(){
        double e = pot->get_energy(x);
        ASSERT_NEAR(e, etrue, 1e-10);
    }

    void test_energy_gradient(){
        g = Array<double>(x.size());
        gnum = Array<double>(g.size());
        double e = pot->get_energy_gradient(x, g);
        EXPECT_NEAR(e, etrue, 1e-10);
        pot->numerical_gradient(x, gnum, 1e-6);
        for (size_t k=0; k<g.size(); ++k){
            EXPECT_NEAR(g[k], gnum[k], 1e-6);
        }
    }

    void test_energy_gradient_hessian(){
        g = Array<double>(x.size());
        gnum = Array<double>(g.size());
        h = Array<double>(x.size()*x.size());
        hnum = Array<double>(h.size());
        double e = pot->get_energy_gradient_hessian(x, g, h);
        double ecomp = pot->get_energy(x);
        pot->numerical_gradient(x, gnum);
        pot->numerical_hessian(x, hnum);

        EXPECT_NEAR(e, ecomp, 1e-10);
        for (size_t i=0; i<g.size(); ++i){
            ASSERT_NEAR(g[i], gnum[i], 1e-6);
        }
        for (size_t i=0; i<h.size(); ++i){
            ASSERT_NEAR(h[i], hnum[i], 1e-3);
        }
    }
};

#endif
