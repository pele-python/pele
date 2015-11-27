#include "pele/array.h"
#include "pele/combination.h"
#include "pele/pspin_spherical.h"
#include "pele/meta_pow.h"
#include "test_utils.hpp"

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <gtest/gtest.h>

static double const EPS = std::numeric_limits<double>::min();
#define EXPECT_NEAR_RELATIVE(A, B, T)  EXPECT_NEAR(A/(fabs(A)+fabs(B) + EPS), B/(fabs(A)+fabs(B) + EPS), T)

using pele::Array;
using pele::MeanFieldPSpinSpherical;
using pele::pos_int_pow;

class MeanFieldPSpinSphericalTest :  public PotentialTest
{
public:
    size_t n;
    Array<double> interactions;
    virtual void setup_potential(){
        pot = std::shared_ptr<pele::BasePotential> (new MeanFieldPSpinSpherical<4>(interactions, n));
    }

    virtual void SetUp(){
        n = 5;
        x = Array<double>(n+1,1);
        interactions = Array<double>(n*n*n*n, 1);
        etrue = -5;
        setup_potential();
    }
};

TEST_F(MeanFieldPSpinSphericalTest, Energy_Works){
    test_energy();
}

TEST_F(MeanFieldPSpinSphericalTest, EnergyGradient_AgreesWithNumerical){
    test_energy_gradient();
}

TEST_F(MeanFieldPSpinSphericalTest, EnergyGradient_AgreesWithNumerical2){
    x[0] = std::sqrt(2);
    x[1] = 3.14159265359;
    x[2] = 1.7;
    x[3] = 1.9;
    x[4] = sqrt(3);
    g = Array<double>(x.size());
    gnum = Array<double>(g.size());
    pot->get_energy_gradient(x, g);
    pot->numerical_gradient(x, gnum, 1e-6);
    for (size_t k=0; k<g.size(); ++k){
        EXPECT_NEAR(g[k], gnum[k], 1e-6);
    }
}

TEST_F(MeanFieldPSpinSphericalTest, EnergyGradientHessian_AgreesWithNumerical){
    test_energy_gradient_hessian();
}

TEST_F(MeanFieldPSpinSphericalTest, EnergyGradientHessian_AgreesWithNumerical2){
    x[0] = std::sqrt(2);
    x[1] = 3.14159265359;
    x[2] = 1.7;
    x[3] = 1.9;
    x[4] = sqrt(3);
    g = Array<double>(x.size());
    gnum = Array<double>(g.size());
    h = Array<double>(x.size()*x.size());
    hnum = Array<double>(h.size());
    double e = pot->get_energy_gradient_hessian(x, g, h);
    double ecomp = pot->get_energy(x);
    pot->numerical_gradient(x, gnum);
    pot->numerical_hessian(x, hnum);

    for (size_t i=0; i<g.size(); ++i){
        ASSERT_NEAR(g[i], gnum[i], 1e-6);
    }
    for (size_t i=0; i<h.size(); ++i){
        ASSERT_NEAR(h[i], hnum[i], 1e-3);
    }
}
