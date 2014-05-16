#include "pele/array.h"
#include "pele/morse.h"
#include "test_utils.hpp"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <gtest/gtest.h>
#include <cmath>

using pele::Array;

class MorseTest :  public PotentialTest
{
public:
    double rho, r0, A;
    void SetUp(){
        rho = 1.2;
        r0 = 1.6;
        A = 1.7;
        x = Array<double>(9);
        x[0]  = 0.1;
        x[1]  = 0.2;
        x[2]  = 0.3;
        x[3]  = 0.44;
        x[4]  = 0.55;
        x[5]  = 1.66;
        x[6] = 0.88;
        x[7] = 1.1;
        x[8] = 3.32;
        etrue = -3.6880170976137334;

        pot = std::shared_ptr<pele::BasePotential> (new pele::Morse (rho, r0, A));
    }
};


TEST_F(MorseTest, Energy_Works){
    test_energy();
}

TEST_F(MorseTest, EnergyGradient_AgreesWithNumerical){
    test_energy_gradient();
}

TEST_F(MorseTest, EnergyGradientHessian_AgreesWithNumerical){
    test_energy_gradient_hessian();
}
