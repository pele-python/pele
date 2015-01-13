#include "pele/array.h"
#include "pele/frozen_atoms.h"
#include "pele/lj.h"

#include <iostream>
#include <stdexcept>
#include <gtest/gtest.h>

using pele::Array;
using pele::LJ;

/**
 * Pairwise Lennard-Jones potential with frozen atoms
 */
class LJFrozen : public pele::FrozenPotentialWrapper{
public:
    LJFrozen(double C6, double C12, Array<double> const & reference_coords, Array<size_t> const & frozen_dof)
        : pele::FrozenPotentialWrapper(std::make_shared<LJ>(C6, C12),
                reference_coords, frozen_dof )
    {}
};


class FrozenLJTest :  public ::testing::Test
{
public:
    double c6, c12, etrue;
    Array<double> x, y;
    Array<size_t> frozen_dof;
    virtual void SetUp(){
        c6 = 1.2;
        c12 = 2.3;
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
        frozen_dof = Array<size_t>(3);
        frozen_dof[0] = 0;
        frozen_dof[1] = 3;
        frozen_dof[2] = 4;
    }
};

TEST_F(FrozenLJTest, TestEnergy_Correct){
    LJ pot_nofreeze(c6, c12);
    LJFrozen pot(c6, c12, x, frozen_dof);
    auto xred = pot.coords_converter.get_reduced_coords(x);
    double e = pot.get_energy(xred);
    double etrue =  pot_nofreeze.get_energy(x);
    EXPECT_NEAR(e, etrue, 1e-10);
}

TEST_F(FrozenLJTest, TestEnergyGradient_Correct){
    LJ pot_nofreeze(c6, c12);
    LJFrozen pot(c6, c12, x, frozen_dof);
    auto xred = pot.coords_converter.get_reduced_coords(x);
    Array<double> gred(xred.size()), gtrue(x.size());
    double e = pot.get_energy_gradient(xred, gred);
    double etrue =  pot_nofreeze.get_energy_gradient(x, gtrue);
    EXPECT_NEAR(e, etrue, 1e-10);
    auto gtrue_red = pot.coords_converter.get_reduced_coords(gtrue);
    for (size_t i=0; i<gred.size(); ++i){
        EXPECT_NEAR(gtrue_red[i], gred[i], 1e-10);
    }
}

TEST_F(FrozenLJTest, TestNumericalGradient_Correct){
    LJ pot_nofreeze(c6, c12);
    LJFrozen pot(c6, c12, x, frozen_dof);
    auto xred = pot.get_reduced_coords(x);
    Array<double> gred(xred.size()), gtrue(x.size());
    pot.numerical_gradient(xred, gred);
    pot_nofreeze.numerical_gradient(x, gtrue);
    auto gtrue_red = pot.get_reduced_coords(gtrue);
    for (size_t i=0; i<gred.size(); ++i){
        EXPECT_NEAR(gtrue_red[i], gred[i], 1e-10);
    }
}

TEST_F(FrozenLJTest, TestNumericalHessian_Correct){
    LJ pot_nofreeze(c6, c12);
    LJFrozen pot(c6, c12, x, frozen_dof);
    auto xred = pot.get_reduced_coords(x);
    Array<double> hred(xred.size()*xred.size()), htrue(x.size()*x.size());
    pot.numerical_hessian(xred, hred);
    pot_nofreeze.numerical_hessian(x, htrue);
    auto htrue_red = pot.coords_converter.get_reduced_hessian(htrue);
//    std::cerr << htrue << "\n";
//    std::cerr << hred << "\n";
//    std::cerr << htrue_red << "\n";
    for (size_t i=0; i<hred.size(); ++i){
        EXPECT_NEAR(htrue_red[i], hred[i], 1e-10);
    }
}

TEST_F(FrozenLJTest, TestEnergyGradientHessian_Correct){
    LJ pot_nofreeze(c6, c12);
    LJFrozen pot(c6, c12, x, frozen_dof);
    auto xred = pot.get_reduced_coords(x);
    Array<double> gred(xred.size()), gtrue(x.size());
    Array<double> hred(xred.size()*xred.size()), htrue(x.size()*x.size());
    double e = pot.get_energy_gradient_hessian(xred, gred, hred);
    double etrue =  pot_nofreeze.get_energy_gradient_hessian(x, gtrue, htrue);
    EXPECT_NEAR(e, etrue, 1e-10);
    auto gtrue_red = pot.get_reduced_coords(gtrue);
    for (size_t i=0; i<gred.size(); ++i){
        EXPECT_NEAR(gtrue_red[i], gred[i], 1e-10);
    }
//    std::cerr << gred << "\n";
//    std::cerr << gtrue_red << "\n";
//    std::cerr << gtrue << "\n";

    auto htrue_red = pot.coords_converter.get_reduced_hessian(htrue);
//    std::cerr << htrue << "\n";
//    std::cerr << hred << "\n";
//    std::cerr << htrue_red << "\n";
    for (size_t i=0; i<hred.size(); ++i){
        EXPECT_NEAR(htrue_red[i], hred[i], 1e-10);
    }
}

TEST_F(FrozenLJTest, TestAccessors){
    LJFrozen pot(c6, c12, x, frozen_dof);
    auto fdof = pot.get_frozen_dof();
    ASSERT_EQ(fdof.size(), frozen_dof.size());
    for (size_t i = 0; i < fdof.size(); ++i) {
        fdof[i] = frozen_dof[i];
    }

}

