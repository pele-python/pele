#include "pele/array.h"
#include "pele/morse.h"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <gtest/gtest.h>
#include <cmath>

using pele::Array;

class MorseTest :  public ::testing::Test
{
public:
    double rho, r0, A, etrue;
    Array<double> x, y;
    void SetUp(){
        rho = 1.2;
        r0 = 1.6;
        A = 1.7;
        x = Array<double>(6);
        x[0]  = 0.1;
        x[1]  = 0.2;
        x[2]  = 0.3;
        x[3]  = 0.44;
        x[4]  = 0.55;
        x[5]  = 1.66;
        etrue = -1.6288465928749187;
    }
};

TEST_F(MorseTest, Energy_Works){
    pele::Morse pot(rho, r0, A);
    double e = pot.get_energy(x);
    ASSERT_NEAR(e, etrue, 1e-10);
}

TEST_F(MorseTest, EnergyGradient_Works){
    pele::Morse pot(rho, r0, A);
    Array<double> g(6);
    Array<double> gnum(6);
    double e = pot.get_energy_gradient(x, g);
    ASSERT_NEAR(e, etrue, 1e-10);
    pot.numerical_gradient(x, gnum, 1e-6);
    for (int k=0; k<6; ++k){
        ASSERT_NEAR(g[k], gnum[k], 1e-6);
    }
}

TEST_F(MorseTest, Hessian_Works){
    pele::Morse pot(rho, r0, A);
    Array<double> h(6*6);
    Array<double> h_num(6*6);
    pot.get_hessian(x, h);
    pot.numerical_hessian(x, h_num);
    for (int i; i<h.size();++i)
        ASSERT_NEAR(h[i], h_num[i],1e-6);
}

//TEST_F(MorseTest, Hessian_Works2){
//    pele::Morse lj(c6, c12);
//    Array<double> h(9*9);
//    Array<double> h_num(9*9);
//    lj.get_hessian(y, h);
//    lj.numerical_hessian(y,h_num);
//    for (int i; i<h.size();++i)
//        ASSERT_NEAR(h[i],h_num[i],1e-6);
//}
