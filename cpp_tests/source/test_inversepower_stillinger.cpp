#include <gtest/gtest.h>

#include "pele/inversepower_stillinger.h"

#include "test_utils.hpp"

static double const EPS = std::numeric_limits<double>::min();
#define EXPECT_NEAR_RELATIVE(A, B, T)  EXPECT_NEAR(A/(fabs(A)+fabs(B) + EPS), B/(fabs(A)+fabs(B) + EPS), T)

class BasicIPSTest : public ::testing::Test {
public:
    size_t ndim;
    size_t ndof;
    size_t npart;
    pele::Array<double> x, radii;
    double etrue;
    size_t exponent;
    std::shared_ptr<pele::InversePowerStillinger<2> > pot;
    double a;
    virtual void SetUp()
    {
        ndim = 2;
        npart = 2;
        ndof = ndim * npart;
        x = {1.1, 3, 1, 2};
        radii = {1.05, 1.05};
        a = radii[0]*2.;
        exponent = 4;
        etrue = std::pow(a / std::sqrt(std::pow(x[0] - x[2], 2) + std::pow(x[1] - x[3], 2)), exponent);
        pot = std::make_shared<pele::InversePowerStillinger<2> >(exponent, radii);
    }
};

TEST_F(BasicIPSTest, Pow_Works)
{
    for (size_t power = 0; power < 129; ++power) {
        const double inp = 0.3;
        pele::InversePowerStillinger_interaction in(power, radii);
        EXPECT_NEAR_RELATIVE(std::pow(inp, power), in.power(inp, power), 1e-14);
    }
    
    const double e = pot->get_energy(x);
    ASSERT_NEAR(etrue, e, 1e-10);
}

class TestInversePowerStillingerAuto : public PotentialTest {
    size_t ndim;
    size_t ndof;
    size_t npart;
    size_t exponent;
    double a;
    pele::Array<double> radii;
    virtual void SetUp()
    {
        ndim = 2;
        npart = 2;
        ndof = ndim * npart;
        x = {1.1, 3, 1, 2};
        radii = {1.05, 1.05};
        a = radii[0]*2;
        exponent = 4;
        etrue = std::pow(a / std::sqrt(std::pow(x[0] - x[2], 2) + std::pow(x[1] - x[3], 2)), exponent);
        pot = std::make_shared<pele::InversePowerStillinger<2> >(exponent, radii);
    }
};

TEST_F(TestInversePowerStillingerAuto, Energy_Works)
{
    test_energy();
}

TEST_F(TestInversePowerStillingerAuto, EnergyGradient_AgreesWithNumerical)
{
    test_energy_gradient();
}

TEST_F(TestInversePowerStillingerAuto, EnergyGradientHessian_AgreesWithNumerical)
{
    test_energy_gradient_hessian();
}
