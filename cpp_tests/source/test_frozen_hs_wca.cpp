#include "pele/array.h"
#include "pele/frozen_atoms.h"
#include "pele/hs_wca.h"
#include "pele/modified_fire.h"

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <gtest/gtest.h>

using pele::Array;
using pele::HS_WCAFrozen;
using pele::HS_WCAPeriodicFrozen;
using pele::HS_WCA;
using pele::HS_WCAPeriodic;
using pele::HS_WCA2D;
using pele::HS_WCA2DFrozen;
using pele::HS_WCAPeriodic2D;
using pele::HS_WCAPeriodic2DFrozen;

#define EXPECT_NEAR_RELATIVE(A, B, T)  EXPECT_NEAR(fabs(A)/(fabs(A)+fabs(B)+1), fabs(B)/(fabs(A)+fabs(B)+1), T)

class FrozenHS_WCATest: public ::testing::Test{
public:
    typedef pele::MODIFIED_FIRE opt_t;
    double eps, sca, etrue;
    Array<double> radii, radii_small, radii_large;
    Array<double> radii2d;
    Array<double> x, y;
    Array<double> x2d;
    Array<size_t> frozen_dof;
    Array<size_t> frozen_dof_2d;
    double* boxvec;
    double* boxvec2d;
    virtual void SetUp(){
        eps=1.0;
        sca=1.2;
        radii.resize(3);
        double f = 0.3; // for too much overlap, numerical accuracy decreases
        radii[0] = 1.0*f;
        radii[1] = 1.1*f;
        radii[2] = 0.9*f;
        radii2d = radii;
        radii_small.resize(3);
        radii_large.resize(3);
        for (size_t i(0); i < radii.size(); ++i){
            radii_small[i] = radii[i]/3.;
            radii_large[i] = radii[i]*1.3;
        }
        boxvec = new double[3];
        boxvec2d = new double[2];
        boxvec[0] = 5;
        boxvec[1] = 6;
        boxvec[2] = 7;
        boxvec2d[0] = 5;
        boxvec2d[1] = 6;
        x.resize(9);
        x[0] = 0.1;
        x[1] = 0.2;
        x[2] = 0.3;
        x[3] = 0.44;
        x[4] = 0.55;
        x[5] = 1.66;
        x[6] = 0.88;
        x[7] = 1.1;
        x[8] = 3.32;
        x2d.resize(6);
        x2d[0] = 0.1;
        x2d[1] = 0.11;
        x2d[2] = 0.6;
        x2d[3] = 0.66;
        x2d[4] = 2.2;
        x2d[5] = 2.3;
        frozen_dof.resize(3);
        frozen_dof[0] = 0;
        frozen_dof[1] = 3;
        frozen_dof[2] = 4;
        frozen_dof_2d = frozen_dof;
    }
    virtual void TearDown() {
        delete[] boxvec;
        delete[] boxvec2d;
    }
};

TEST_F(FrozenHS_WCATest, TestEnergy_Correct){
    HS_WCA pot_nofreeze(eps,sca,radii);
    HS_WCAFrozen pot(eps, sca, radii, x, frozen_dof);
    auto xred = pot.coords_converter.get_reduced_coords(x);
    double e = pot.get_energy(xred);
    double etrue =  pot_nofreeze.get_energy(x);
    //std::cout << "etrue: " << etrue << std::endl;
    EXPECT_NEAR(e, etrue, 1e-10);
}

TEST_F(FrozenHS_WCATest, TestEnergy_Correct_2D){
    HS_WCA2D pot_nofreeze(eps,sca,radii2d);
    HS_WCA2DFrozen pot(eps, sca, radii2d, x2d, frozen_dof_2d);
    auto xred = pot.coords_converter.get_reduced_coords(x2d);
    double e = pot.get_energy(xred);
    double etrue =  pot_nofreeze.get_energy(x2d);
    //std::cout << "etrue: " << etrue << std::endl;
    EXPECT_NEAR(e, etrue, 1e-10);
}

TEST_F(FrozenHS_WCATest, TestEnergy_Correct_Periodic){
    HS_WCAPeriodic pot_nofreeze(eps, sca, radii, boxvec);
    HS_WCAPeriodicFrozen pot(eps, sca, radii, boxvec, x, frozen_dof);
    auto xred = pot.coords_converter.get_reduced_coords(x);
    double e = pot.get_energy(xred);
    double etrue =  pot_nofreeze.get_energy(x);
    EXPECT_NEAR(e, etrue, 1e-10);
}

TEST_F(FrozenHS_WCATest, TestRepulsive_Correct){
    HS_WCA pot_nofreeze(eps,sca,radii);
    HS_WCAFrozen pot(eps, sca, radii, x, frozen_dof);
    HS_WCAFrozen pot_small(eps, sca, radii_small, x, frozen_dof);
    HS_WCAFrozen pot_large(eps, sca, radii_large, x, frozen_dof);
    auto xred = pot.coords_converter.get_reduced_coords(x);
    double e = pot.get_energy(xred);
    double e_small = pot_small.get_energy(xred);
    double e_large = pot_large.get_energy(xred);
    double etrue =  pot_nofreeze.get_energy(x);
    EXPECT_NEAR(e, etrue, 1e-10);
    EXPECT_TRUE(e_small<=e);
    EXPECT_TRUE(e<=e_large);
}

TEST_F(FrozenHS_WCATest, TestRepulsive_Correct_Periodic){
    HS_WCAPeriodic pot_nofreeze(eps, sca, radii, boxvec);
    HS_WCAPeriodicFrozen pot(eps, sca, radii, boxvec, x, frozen_dof);
    HS_WCAPeriodicFrozen pot_small(eps, sca, radii_small, boxvec, x, frozen_dof);
    HS_WCAPeriodicFrozen pot_large(eps, sca, radii_large, boxvec, x, frozen_dof);
    auto xred = pot.coords_converter.get_reduced_coords(x);
    double e = pot.get_energy(xred);
    double e_small = pot_small.get_energy(xred);
    double e_large = pot_large.get_energy(xred);
    double etrue =  pot_nofreeze.get_energy(x);
    EXPECT_NEAR(e, etrue, 1e-10);
    EXPECT_TRUE(e_small<=e);
    EXPECT_TRUE(e<=e_large);
}

TEST_F(FrozenHS_WCATest, TestEnergyGradient_Correct){
    HS_WCA pot_nofreeze(eps,sca,radii);
    HS_WCAFrozen pot(eps, sca, radii, x, frozen_dof);
    auto xred = pot.coords_converter.get_reduced_coords(x);
    Array<double> gred(xred.size()), gtrue(x.size());
    double e = pot.get_energy_gradient(xred,gred);
    double etrue = pot_nofreeze.get_energy_gradient(x,gtrue);
    EXPECT_NEAR(e, etrue, 1e-10);
    auto gtrue_red = pot.coords_converter.get_reduced_coords(gtrue);
    for (size_t i=0; i<gred.size(); ++i){
        EXPECT_NEAR(gtrue_red[i], gred[i], 1e-10);
    }
}

TEST_F(FrozenHS_WCATest, TestEnergyGradient_Correct_2D){
    HS_WCA2D pot_nofreeze(eps,sca,radii2d);
    HS_WCA2DFrozen pot(eps, sca, radii2d, x2d, frozen_dof_2d);
    auto xred = pot.coords_converter.get_reduced_coords(x2d);
    Array<double> gred(xred.size()), gtrue(x2d.size());
    double e = pot.get_energy_gradient(xred,gred);
    double etrue = pot_nofreeze.get_energy_gradient(x2d,gtrue);
    EXPECT_NEAR(e, etrue, 1e-10);
    auto gtrue_red = pot.coords_converter.get_reduced_coords(gtrue);
    for (size_t i=0; i<gred.size(); ++i){
        EXPECT_NEAR(gtrue_red[i], gred[i], 1e-10);
    }
}

TEST_F(FrozenHS_WCATest, TestEnergyGradient_Correct_Periodic){
    HS_WCAPeriodic pot_nofreeze(eps, sca, radii, boxvec);
    HS_WCAPeriodicFrozen pot(eps, sca, radii, boxvec, x, frozen_dof);
    auto xred = pot.coords_converter.get_reduced_coords(x);
    Array<double> gred(xred.size()), gtrue(x.size());
    double e = pot.get_energy_gradient(xred,gred);
    double etrue = pot_nofreeze.get_energy_gradient(x,gtrue);
    EXPECT_NEAR(e, etrue, 1e-10);
    auto gtrue_red = pot.coords_converter.get_reduced_coords(gtrue);
    for (size_t i=0; i<gred.size(); ++i){
        EXPECT_NEAR(gtrue_red[i], gred[i], 1e-10);
    }
}

TEST_F(FrozenHS_WCATest, TestNumericalGradient_Correct){
    HS_WCA pot_nofreeze(eps,sca,radii);
    HS_WCAFrozen pot(eps, sca, radii, x, frozen_dof);
    auto xred = pot.coords_converter.get_reduced_coords(x);
    Array<double> gred(xred.size()), gtrue(x.size());
    pot.numerical_gradient(xred,gred);
    pot_nofreeze.numerical_gradient(x,gtrue);
    auto gtrue_red = pot.coords_converter.get_reduced_coords(gtrue);
    for (size_t i=0; i<gred.size(); ++i){
        EXPECT_NEAR(gtrue_red[i], gred[i], 1e-10);
    }
}

TEST_F(FrozenHS_WCATest, TestNumericalGradient_Correct_Periodic){
    HS_WCAPeriodic pot_nofreeze(eps, sca, radii, boxvec);
    HS_WCAPeriodicFrozen pot(eps, sca, radii, boxvec, x, frozen_dof);
    auto xred = pot.coords_converter.get_reduced_coords(x);
    Array<double> gred(xred.size()), gtrue(x.size());
    pot.numerical_gradient(xred,gred);
    pot_nofreeze.numerical_gradient(x,gtrue);
    auto gtrue_red = pot.coords_converter.get_reduced_coords(gtrue);
    for (size_t i=0; i<gred.size(); ++i){
        EXPECT_NEAR(gtrue_red[i], gred[i], 1e-10);
    }
}

TEST_F(FrozenHS_WCATest, TestNumericalHessian_Correct){
    HS_WCA pot_nofreeze(eps, sca, radii);
    HS_WCAFrozen pot(eps, sca, radii, x, frozen_dof);
    auto xred = pot.coords_converter.get_reduced_coords(x);
    Array<double> hred(xred.size()*xred.size()), htrue(x.size()*x.size());
    pot.numerical_hessian(xred, hred);
    pot_nofreeze.numerical_hessian(x, htrue);
    auto htrue_red = pot.coords_converter.get_reduced_hessian(htrue);
    for (size_t i=0; i<hred.size(); ++i){
        EXPECT_NEAR(htrue_red[i], hred[i], 1e-10);
    }
}

TEST_F(FrozenHS_WCATest, TestNumericalHessian_Correct_Periodic){
    HS_WCAPeriodic pot_nofreeze(eps, sca, radii, boxvec);
    HS_WCAPeriodicFrozen pot(eps, sca, radii, boxvec, x, frozen_dof);
    auto xred = pot.coords_converter.get_reduced_coords(x);
    Array<double> hred(xred.size()*xred.size()), htrue(x.size()*x.size());
    pot.numerical_hessian(xred, hred);
    pot_nofreeze.numerical_hessian(x, htrue);
    auto htrue_red = pot.coords_converter.get_reduced_hessian(htrue);
    for (size_t i=0; i<hred.size(); ++i){
        EXPECT_NEAR(htrue_red[i], hred[i], 1e-10);
    }
}

TEST_F(FrozenHS_WCATest, TestEnergyGradientHessian_Correct){
    HS_WCA pot_nofreeze(eps, sca, radii);
    HS_WCAFrozen pot(eps, sca, radii, x, frozen_dof);
    //LJ pot_nofreeze(1.2, 2.3);
    //LJFrozen pot(1.2, 2.3, x, frozen_dof);
    auto xred = pot.coords_converter.get_reduced_coords(x);
    Array<double> gred(xred.size()), gtrue(x.size());
    Array<double> hred(xred.size()*xred.size()), htrue(x.size()*x.size());
    double e = pot.get_energy_gradient_hessian(xred, gred, hred);
    double etrue =  pot_nofreeze.get_energy_gradient_hessian(x, gtrue, htrue);
    EXPECT_NEAR(e, etrue, 1e-10);
    auto gtrue_red = pot.coords_converter.get_reduced_coords(gtrue);
    for (size_t i=0; i<gred.size(); ++i){
        EXPECT_NEAR(gtrue_red[i], gred[i], 1e-10);
    }
    auto htrue_red = pot.coords_converter.get_reduced_hessian(htrue);
    for (size_t i=0; i<hred.size(); ++i){
        EXPECT_NEAR(htrue_red[i], hred[i], 1e-10);
    }
}

TEST_F(FrozenHS_WCATest, TestEnergyGradientHessian_Correct_2D){
    HS_WCA2D pot_nofreeze(eps, sca, radii2d);
    HS_WCA2DFrozen pot(eps, sca, radii2d, x2d, frozen_dof_2d);
    auto xred = pot.coords_converter.get_reduced_coords(x2d);
    Array<double> gred(xred.size()), gtrue(x2d.size());
    Array<double> hred(xred.size()*xred.size()), htrue(x2d.size()*x2d.size());
    double e = pot.get_energy_gradient_hessian(xred, gred, hred);
    double etrue =  pot_nofreeze.get_energy_gradient_hessian(x2d, gtrue, htrue);
    EXPECT_NEAR(e, etrue, 1e-10);
    auto gtrue_red = pot.coords_converter.get_reduced_coords(gtrue);
    for (size_t i=0; i<gred.size(); ++i){
        EXPECT_NEAR(gtrue_red[i], gred[i], 1e-10);
    }
    auto htrue_red = pot.coords_converter.get_reduced_hessian(htrue);
    for (size_t i=0; i<hred.size(); ++i){
        EXPECT_NEAR(htrue_red[i], hred[i], 1e-10);
    }
}

TEST_F(FrozenHS_WCATest, TestEnergyGradientHessian_Correct_Periodic){
    HS_WCAPeriodic pot_nofreeze(eps, sca, radii, boxvec);
    HS_WCAPeriodicFrozen pot(eps, sca, radii, boxvec, x, frozen_dof);
    auto xred = pot.coords_converter.get_reduced_coords(x);
    Array<double> gred(xred.size()), gtrue(x.size());
    Array<double> hred(xred.size()*xred.size()), htrue(x.size()*x.size());
    double e = pot.get_energy_gradient_hessian(xred, gred, hred);
    double etrue =  pot_nofreeze.get_energy_gradient_hessian(x, gtrue, htrue);
    EXPECT_NEAR(e, etrue, 1e-10);
    auto gtrue_red = pot.coords_converter.get_reduced_coords(gtrue);
    for (size_t i=0; i<gred.size(); ++i){
        EXPECT_NEAR(gtrue_red[i], gred[i], 1e-10);
    }
    auto htrue_red = pot.coords_converter.get_reduced_hessian(htrue);
    for (size_t i=0; i<hred.size(); ++i){
        EXPECT_NEAR(htrue_red[i], hred[i], 1e-10);
    }
}

TEST_F(FrozenHS_WCATest, TestEnergyGradientHessian_AgreesWithNumerical){
    HS_WCAFrozen pot(eps, sca, radii, x, frozen_dof);
    auto xred = pot.coords_converter.get_reduced_coords(x);
    Array<double> grad(xred.size());
    Array<double> grad_num(grad.size());
    Array<double> hess(xred.size()*xred.size());
    Array<double> hess_num(hess.size());
    const double energy = pot.get_energy_gradient_hessian(xred, grad, hess);
    const double energy_comp(pot.get_energy(xred));
    pot.numerical_gradient(xred, grad_num);
    pot.numerical_hessian(xred, hess_num);
    EXPECT_NEAR(energy, energy_comp, 1e-10);
    for (size_t i(0); i < grad.size(); ++i){
        EXPECT_NEAR(grad[i], grad_num[i], 1e-10);
    }
    for (size_t i(0); i < hess.size(); ++i){
        EXPECT_NEAR(hess[i], hess_num[i], 1e-5);
    }
}

TEST_F(FrozenHS_WCATest, TestEnergyGradientHessian_AgreesWithNumerical_2D){
    HS_WCA2DFrozen pot(eps, sca, radii2d, x2d, frozen_dof_2d);
    auto xred = pot.coords_converter.get_reduced_coords(x2d);
    Array<double> grad(xred.size());
    Array<double> grad_num(grad.size());
    Array<double> hess(xred.size()*xred.size());
    Array<double> hess_num(hess.size());
    const double energy = pot.get_energy_gradient_hessian(xred, grad, hess);
    const double energy_comp(pot.get_energy(xred));
    pot.numerical_gradient(xred, grad_num);
    pot.numerical_hessian(xred, hess_num);
    EXPECT_NEAR(energy, energy_comp, 1e-10);
    for (size_t i(0); i < grad.size(); ++i){
	EXPECT_NEAR_RELATIVE(grad[i], grad_num[i], 1e-8);
    }
    for (size_t i(0); i < hess.size(); ++i){
	EXPECT_NEAR_RELATIVE(hess[i], hess_num[i], 1e-8);
    }
}

TEST_F(FrozenHS_WCATest, TestMinimizationFreezing_Correct){
    // minimization for non-frozen HS_WCA
    HS_WCA pot_nofreeze(eps,sca,radii);
    const auto e_notfreeze_before = pot_nofreeze.get_energy(x);
    opt_t minimizer(&pot_nofreeze, x, 1e-4, 1.0, 0.5);
    minimizer.run();
    const auto x_after = minimizer.get_x();
    const auto e_notfreeze_after = pot_nofreeze.get_energy(x_after);
    EXPECT_TRUE( e_notfreeze_after <= e_notfreeze_before );
    EXPECT_NEAR( e_notfreeze_after, 0, 1e-10 );
    // minimization for frozen HS_WCA
    HS_WCAFrozen pot(eps, sca, radii, x, frozen_dof);
    auto xred = pot.coords_converter.get_reduced_coords(x);
    const auto e_before = pot.get_energy(xred);
    opt_t minimizer_red(&pot, xred, 1e-4, 1.0, 0.5);
    minimizer_red.run();
    const auto xred_after = minimizer_red.get_x();
    const auto e_after = pot.get_energy(xred_after);
    EXPECT_TRUE( e_after <= e_before );
    EXPECT_TRUE( e_notfreeze_after <= e_after );
    // verify that frozen dof are in fact frozen
    const auto xred_after_inflated = pot.coords_converter.get_full_coords(xred_after);
    for (size_t i(0); i < frozen_dof.size(); ++i){
        const auto idx = frozen_dof[i];
        EXPECT_NEAR( x[idx], xred_after_inflated[idx], 1e-10 );
    }
}

TEST_F(FrozenHS_WCATest, TestMinimizationFreezing_Correct_2D){
    HS_WCA2D pot_nofreeze(eps,sca,radii2d);
    const auto e_notfreeze_before = pot_nofreeze.get_energy(x2d);
    opt_t minimizer(&pot_nofreeze, x2d, 1e-4, 1.0, 0.5);
    minimizer.run();
    const auto x_after = minimizer.get_x();
    const auto e_notfreeze_after = pot_nofreeze.get_energy(x_after);
    EXPECT_TRUE( e_notfreeze_after <= e_notfreeze_before );
    EXPECT_NEAR( e_notfreeze_after, 0, 1e-10 );
    HS_WCA2DFrozen pot(eps, sca, radii2d, x2d, frozen_dof_2d);
    auto xred = pot.coords_converter.get_reduced_coords(x2d);
    const auto e_before = pot.get_energy(xred);
    opt_t minimizer_red(&pot, xred, 1e-4, 1.0, 0.5);
    minimizer_red.run();
    const auto xred_after = minimizer_red.get_x();
    const auto e_after = pot.get_energy(xred_after);
    EXPECT_TRUE( e_after <= e_before );
    EXPECT_TRUE( e_notfreeze_after <= e_after );
    // verify that frozen dof are in fact frozen
    const auto xred_after_inflated = pot.coords_converter.get_full_coords(xred_after);
    for (size_t i(0); i < frozen_dof_2d.size(); ++i){
	const auto idx = frozen_dof_2d[i];
	EXPECT_NEAR( x2d[idx], xred_after_inflated[idx], 1e-10 );
    }
}

