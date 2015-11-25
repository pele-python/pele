#include "pele/array.h"
#include "pele/combination.h"
#include "pele/mf_p_spin_spherical.h"
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
        n = 20;
        x = Array<double>(n+1,1);
        interactions = Array<double>(n*n*n*n, 1);
        etrue = -4845;
        setup_potential();
    }
};

TEST_F(MeanFieldPSpinSphericalTest, Energy_Works){
    test_energy();
}

TEST_F(MeanFieldPSpinSphericalTest, EnergyGradient_AgreesWithNumerical){
    test_energy_gradient();
}

TEST_F(MeanFieldPSpinSphericalTest, EnergyGradientHessian_AgreesWithNumerical){
    test_energy_gradient_hessian();
}

