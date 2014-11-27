#include "pele/array.h"
#include "pele/wca.h"
#include "pele/hs_wca.h"
#include "pele/neighbor_iterator.h"
#include "test_utils.hpp"

#include <iostream>
#include <stdexcept>
#include <gtest/gtest.h>

using pele::Array;
using pele::WCA;
using pele::HS_WCA;


class WCATest :  public PotentialTest
{
public:
    double sig, eps;
    size_t natoms;

    virtual void setup_potential(){
        pot = std::shared_ptr<pele::BasePotential> (new pele::WCA(sig, eps));
    }

    virtual void SetUp(){
        natoms = 3;
        sig = 1.4;
        eps = 2.1;
        x = Array<double>(natoms*3);
        x[0] = 0.1;
        x[1] = 0.2;
        x[2] = 0.3;
        x[3] = 0.44;
        x[4] = 0.55;
        x[5] = 1.66;
        x[6] = 0.88;
        x[7] = 1.1;
        x[8] = 3.32;
        etrue = 0.9009099166892105;

        setup_potential();
    }
};

TEST_F(WCATest, Energy_Works){
    test_energy();
}
TEST_F(WCATest, EnergyGradient_AgreesWithNumerical){
    test_energy_gradient();
}
TEST_F(WCATest, EnergyGradientHessian_AgreesWithNumerical){
    test_energy_gradient_hessian();
}

class WCAAtomListTest :  public WCATest
{
public:
    virtual void setup_potential(){
        pele::Array<size_t> atoms(natoms);
        for (size_t i =0; i<atoms.size(); ++i){
            atoms[i] = i;
        }
        pot = std::shared_ptr<pele::BasePotential> (new pele::WCAAtomList(
                sig, eps, atoms
                ));
    }
};

TEST_F(WCAAtomListTest, Energy_Works){
    test_energy();
}
TEST_F(WCAAtomListTest, EnergyGradient_AgreesWithNumerical){
    test_energy_gradient();
}
TEST_F(WCAAtomListTest, EnergyGradientHessian_AgreesWithNumerical){
    test_energy_gradient_hessian();
}



/*
 * HS_WCA tests
 */
class HS_WCATest :  public ::testing::Test
{
public:
    double eps, sca, etrue;
    Array<double> x, g, gnum, radii;
    virtual void SetUp(){
        eps = 2.1;
        sca = 1.4;
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
        double f = .35;
        radii[0] = .91 * f;
        radii[1] = 1.1 * f;
        radii[2] = 1.13 * f;
        etrue = 189.41835811474974;
        g = Array<double>(x.size());
        gnum = Array<double>(x.size());
    }
};

TEST_F(HS_WCATest, Energy_Works){
    HS_WCA<3> pot(eps, sca, radii);
    double e = pot.get_energy(x);
    ASSERT_NEAR(e, etrue, 1e-10);
}

TEST_F(HS_WCATest, EnergyGradient_AgreesWithNumerical){
    HS_WCA<3> pot(eps, sca, radii);
    double e = pot.get_energy_gradient(x, g);
    double ecomp = pot.get_energy(x);
    ASSERT_NEAR(e, ecomp, 1e-10);
    pot.numerical_gradient(x, gnum, 1e-6);
    for (size_t k=0; k<6; ++k){
        ASSERT_NEAR(g[k], gnum[k], 1e-6);
    }
}

TEST_F(HS_WCATest, EnergyGradientHessian_AgreesWithNumerical){
    HS_WCA<3> pot(eps, sca, radii);
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
