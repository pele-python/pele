#include <iostream>
#include <stdexcept>
#include <gtest/gtest.h>

#include "pele/array.h"
#include "pele/inversepower.h"

using pele::Array;
using pele::InversePower;

/*
 * InversePower tests
 */

class InversePowerTest :  public ::testing::Test
{
public:
    double pow, eps, etrue;
    Array<double> x, g, gnum, radii;
    virtual void SetUp(){
    	pow = 2.5;
    	eps = 1;
        x = Array<double>(9);
        x[0] = 0.1;
        x[1] = 0.2;
        x[2] = 0.3;
        x[3] = 0.44;
        x[4] = 0.55;
        x[5] = 1.66;
        x[6] = 0.88;
        x[7] = 1.1;
        x[8] = 3.32;
        radii = Array<double>(3);
        double f = 1.;
        radii[0] = .91 * f;
        radii[1] = 1.1 * f;
        radii[2] = 1.13 * f;
        etrue = 0.023173380132354017;
        g = Array<double>(x.size());
        gnum = Array<double>(x.size());
    }
};

TEST_F(InversePowerTest, Energy_Works){
    InversePower<3> pot(pow, eps, radii);
    double e = pot.get_energy(x);
    ASSERT_NEAR(e, etrue, 1e-10);
}

TEST_F(InversePowerTest, EnergyGradient_AgreesWithNumerical){
	InversePower<3> pot(pow, eps, radii);
	double e = pot.get_energy_gradient(x, g);
	std::cout<<"energy"<<e<<std::endl;
    double ecomp = pot.get_energy(x);
    ASSERT_NEAR(e, ecomp, 1e-10);
    pot.numerical_gradient(x, gnum, 1e-6);
    for (size_t k=0; k<6; ++k){
        ASSERT_NEAR(g[k], gnum[k], 1e-6);
    }
}

TEST_F(InversePowerTest, EnergyGradientHessian_AgreesWithNumerical){
	InversePower<3> pot(pow, eps, radii);
    Array<double> h(x.size()*x.size());
    Array<double> hnum(h.size());
    double e = pot.get_energy_gradient_hessian(x, g, h);
    double ecomp = pot.get_energy(x);
    pot.numerical_gradient(x, gnum);
    pot.numerical_hessian(x, hnum);

    EXPECT_NEAR(e, ecomp, 1e-10);
    for (size_t i=0; i<g.size(); ++i){
        ASSERT_NEAR(g[i], gnum[i], 1e-6);
    }
    for (size_t i=0; i<h.size(); ++i){
        ASSERT_NEAR(h[i], hnum[i], 1e-3);
    }
}

TEST_F(InversePowerTest, MetaPowFunctionsBasic_Work){
    const double op = 42.42;
    const int POW= 5;
    double true_result_direct = op * op * op * op * op;
    double true_result_std = std::pow(op, POW);
    EXPECT_DOUBLE_EQ(true_result_direct, true_result_std);
    EXPECT_DOUBLE_EQ(true_result_direct, pele::pos_int_pow<POW>(op));
    true_result_direct = double(1) / true_result_direct;
    true_result_std = std::pow(op, - POW);
    EXPECT_DOUBLE_EQ(true_result_direct, true_result_std);
    EXPECT_DOUBLE_EQ(true_result_direct, pele::neg_int_pow<- POW>(op));
    true_result_direct = std::sqrt(op * op * op * op * op);
    true_result_std = std::pow(op, 0.5 * POW);
    EXPECT_DOUBLE_EQ(true_result_direct, true_result_std);
    EXPECT_DOUBLE_EQ(true_result_direct, pele::pos_half_int_pow<POW>(op));
    true_result_direct = double(1) / std::sqrt(op * op * op * op * op);
    true_result_std = std::pow(op, - 0.5 * POW);
    EXPECT_DOUBLE_EQ(true_result_direct, true_result_std);
    EXPECT_DOUBLE_EQ(true_result_direct, pele::neg_half_int_pow<- POW>(op));
}

/*
TEST_F(InversePowerTest, MetaPowFunctionsLoop_Work){

}
*/
