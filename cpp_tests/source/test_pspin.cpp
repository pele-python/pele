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
using pele::Factorial;

class MeanFieldPSpinSphericalTestP4 :  public PotentialTest
{
public:
    size_t n;
    Array<double> interactions;
    double tol;
    virtual void setup_potential(){
        pot = std::shared_ptr<pele::BasePotential> (new MeanFieldPSpinSpherical<4>(interactions, n, tol));
    }

    virtual void SetUp(){
        tol=1e-10;
        n = 5;
        x = Array<double>(n,1);
        interactions = Array<double>(n*n*n*n, 1);
        etrue = -5/std::pow(n,3./2.);
        setup_potential();
    }
};

class MeanFieldPSpinSphericalTestP4_2 :  public PotentialTest
{
public:
    size_t n;
    Array<double> interactions;
    double tol;
    virtual void setup_potential(){
        pot = std::shared_ptr<pele::BasePotential> (new MeanFieldPSpinSpherical<4>(interactions, n, tol));
    }

    virtual void SetUp(){
        tol=1e-10;
        n = 5;
        x = Array<double>(n);
        x[0] = std::sqrt(2);
        x[1] = 3.14159265359;
        x[2] = 1.7;
        x[3] = 1.9;
        x[4] = sqrt(3);
        x /= (norm(x)/sqrt(n));
        interactions = Array<double>(n*n*n*n, 1);
        etrue = -3.6975625780267314/std::pow(n,3./2.);
        setup_potential();
    }
};

class MeanFieldPSpinSphericalTestP2 :  public PotentialTest
{
public:
    size_t n;
    Array<double> interactions;
    double tol;
    virtual void setup_potential(){
        pot = std::shared_ptr<pele::BasePotential> (new MeanFieldPSpinSpherical<2>(interactions, n, tol));
    }

    virtual void SetUp(){
        tol=1e-10;
        n = 5;
        x = Array<double>(n,1);
        interactions = Array<double>(n*n, 1);
        etrue = -10;
        setup_potential();
    }
};

class MeanFieldPSpinSphericalTestP2_2 :  public PotentialTest
{
public:
    size_t n;
    Array<double> interactions;
    double tol;
    virtual void setup_potential(){
        pot = std::shared_ptr<pele::BasePotential> (new MeanFieldPSpinSpherical<2>(interactions, n, tol));
    }

    virtual void SetUp(){
        tol=1e-10;
        n = 5;
        x = Array<double>(n);
        x[0] = std::sqrt(2);
        x[1] = 3.14159265359;
        x[2] = 1.7;
        x[3] = 1.9;
        x[4] = sqrt(3);
        x /= (norm(x)/sqrt(n));
        interactions = Array<double>(n*n, 1);
        etrue = -8.9379417937213432;
        setup_potential();
    }
};

TEST_F(MeanFieldPSpinSphericalTestP4, Energy_Works){
    test_energy();
}

TEST_F(MeanFieldPSpinSphericalTestP4, EnergyGradient_AgreesWithNumerical){
    test_energy_gradient();
}

TEST_F(MeanFieldPSpinSphericalTestP4, EnergyGradientHessian_AgreesWithNumerical){
    test_energy_gradient_hessian();
}

TEST_F(MeanFieldPSpinSphericalTestP4_2, Energy_Works){
    test_energy();
}

TEST_F(MeanFieldPSpinSphericalTestP4_2, EnergyGradient_AgreesWithNumerical){
    test_energy_gradient();
}

TEST_F(MeanFieldPSpinSphericalTestP4_2, EnergyGradientHessian_AgreesWithNumerical){
    test_energy_gradient_hessian();
}

TEST_F(MeanFieldPSpinSphericalTestP2, Energy_Works){
    test_energy();
}

TEST_F(MeanFieldPSpinSphericalTestP2, EnergyGradient_AgreesWithNumerical){
    test_energy_gradient();
}

TEST_F(MeanFieldPSpinSphericalTestP2, EnergyGradientHessian_AgreesWithNumerical){
    test_energy_gradient_hessian();
}

TEST_F(MeanFieldPSpinSphericalTestP2_2, Energy_Works){
    test_energy();
}

TEST_F(MeanFieldPSpinSphericalTestP2_2, EnergyGradient_AgreesWithNumerical){
    test_energy_gradient();
}

TEST_F(MeanFieldPSpinSphericalTestP2_2, EnergyGradientHessian_AgreesWithNumerical){
    test_energy_gradient_hessian();
}
