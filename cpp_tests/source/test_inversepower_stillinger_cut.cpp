#include <gtest/gtest.h>

#include "pele/inversepower_stillinger_cut.h"

#include "test_utils.hpp"

static double const EPS = std::numeric_limits<double>::min();
#define EXPECT_NEAR_RELATIVE(A, B, T)  EXPECT_NEAR(A/(fabs(A)+fabs(B) + EPS), B/(fabs(A)+fabs(B) + EPS), T)

class BasicIPSCutTest : public ::testing::Test {
public:
    size_t ndim;
    size_t ndof;
    size_t npart;
    pele::Array<double> x, radii;
    size_t exponent;
    std::shared_ptr<pele::InversePowerStillingerCut<2> > pot;
    double a, rcut, eps;
    double etrue, e0;
    virtual void SetUp()
    {
        ndim = 2;
        npart = 2;
        rcut = 1.5;
        eps = 1e-15;
        ndof = ndim * npart;
        x = {1.1, 3, 1, 2};
        radii = {1.05, 1.05};
        a = radii[0]*2.;
        exponent = 4;
        double dr = std::sqrt(std::pow(x[0] - x[2], 2) + std::pow(x[1] - x[3], 2));
        etrue = get_test_energy(dr, rcut);
        std::cout<<"test energy "<<etrue<<std::endl;
        e0 = get_test_energy(rcut-eps, rcut);
        pot = std::make_shared<pele::InversePowerStillingerCut<2> >(exponent, radii, rcut);
    }

    double get_test_energy(double dr, double rcut){
        double etest = std::pow(a / dr, exponent) + std::pow(a/rcut, exponent) * (
        - 0.5 * (exponent+2) * (exponent+1)
        + dr * exponent * (exponent+2)/rcut
        - 0.5 * dr * dr * exponent*(exponent+1)/(rcut*rcut)
        );
        return etest;
    }
};

TEST_F(BasicIPSCutTest, Pow_Works)
{
    for (size_t power = 0; power < 129; ++power) {
        const double inp = 0.3;
        pele::InversePowerStillinger_cut_interaction in(power, radii, rcut);
        EXPECT_NEAR_RELATIVE(std::pow(inp, power), in.power(inp, power), 1e-14);
    }
    
    const double e = pot->get_energy(x);
    ASSERT_NEAR(etrue, e, 1e-10);
}

TEST_F(BasicIPSCutTest, e0_works)
{
    double dr0 = rcut - eps;
    pele::Array<double> x0 = {0, 0, dr0, 0};
    const double e = pot->get_energy(x0);
    ASSERT_NEAR(e, e0, 1e-14);
    ASSERT_NEAR(e, 0, 1e-14);
}

class TestInversePowerStillingerCutAuto : public PotentialTest {
    size_t ndim;
    size_t ndof;
    size_t npart;
    size_t exponent;
    double a, rcut;
    //do not declare etrue here
    pele::Array<double> radii;
    virtual void SetUp()
    {
        ndim = 2;
        npart = 2;
        ndof = ndim * npart;
        rcut = 1.5;
        x = {1.1, 3, 1, 2};
        radii = {1.05, 1.05};
        a = radii[0]*2;
        exponent = 4;
        double dr = std::sqrt(std::pow(x[0] - x[2], 2) + std::pow(x[1] - x[3], 2));
        etrue = get_test_energy(dr, rcut);
        std::cout<<"test energy "<<etrue<<std::endl;
        pot = std::make_shared<pele::InversePowerStillingerCut<2> >(exponent, radii, rcut);
    }

    double get_test_energy(double dr, double rcut){
        double etest = std::pow(a / dr, exponent) + std::pow(a/rcut, exponent) * (
        - 0.5 * (exponent+2) * (exponent+1)
        + dr * exponent * (exponent+2)/rcut
        - 0.5 * dr * dr * exponent*(exponent+1)/(rcut*rcut)
        );
        return etest;
    }
};

TEST_F(TestInversePowerStillingerCutAuto, Energy_Works)
{
    test_energy();
}

TEST_F(TestInversePowerStillingerCutAuto, EnergyGradient_AgreesWithNumerical)
{
    test_energy_gradient();
}

TEST_F(TestInversePowerStillingerCutAuto, EnergyGradientHessian_AgreesWithNumerical)
{
    test_energy_gradient_hessian();
}
