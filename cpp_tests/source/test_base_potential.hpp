#include "pele/array.h"
#include "pele/base_potential.h"

#include <iostream>
#include <stdexcept>
#include <vector>
#include <gtest/gtest.h>
#include <cmath>

using pele::Array;
using pele::BasePotential;

class HarmonicE : public pele::BasePotential
{
public:
    double get_energy(Array<double> x){
        double energy = 0;
        for (int k=0; k<x.size(); ++k){
            energy += x[k] * x[k];
        }
        return energy / 2.;
    }
};

class BasePotentialTest :  public ::testing::Test
{
public:
    double etrue;
    Array<double> x, g, gtrue, hess, htrue;
    virtual void SetUp(){
        x.resize(2);
        g.resize(2);
        hess.resize(2*2);
        htrue = Array<double>(4, 0.);
        x[0] = 1.;
        x[1] = 1.5;
        etrue = (x[0] * x[0] + x[1] * x[1]) / 2.;
        gtrue = x.copy();
        htrue[0] = 1.;
        htrue[3] = 1.;
    }
};


TEST_F(BasePotentialTest, EOnlyEnergy_Works){
    HarmonicE pot;
    BasePotential * ptr = &pot;
    double e = pot.get_energy(x);
    EXPECT_NEAR(e, etrue, 1e-10);
    double e_ptr = ptr->get_energy(x);
    EXPECT_NEAR(e_ptr, etrue, 1e-10);
}

TEST_F(BasePotentialTest, EOnlyGrad_Works){
    HarmonicE pot;
    BasePotential * ptr = &pot;
    double e = pot.get_energy_gradient(x, g);
    EXPECT_NEAR(e, etrue, 1e-10);

    // the gradient is computed numerically
    for (int k=0; k<x.size(); ++k){
        EXPECT_NEAR(g[k], gtrue[k], 1e-6);
    }
}

TEST_F(BasePotentialTest, EOnlyHess_Works){
    HarmonicE pot;
    BasePotential * ptr = &pot;
    double e = pot.get_energy_gradient_hessian(x, g, hess);
    EXPECT_NEAR(e, etrue, 1e-10);

    // the gradient is computed numerically
    for (int k=0; k<x.size(); ++k){
        EXPECT_NEAR(g[k], gtrue[k], 1e-6);
    }

    // the hessian is computed numerically
    for (int k=0; k<hess.size(); ++k){
        EXPECT_NEAR(hess[k], htrue[k], 1e-3);
    }
}

TEST_F(BasePotentialTest, EOnlyGetHess_Works){
    HarmonicE pot;
    BasePotential * ptr = &pot;
    pot.get_hessian(x, hess);

    // the hessian is computed numerically
    for (int k=0; k<hess.size(); ++k){
        EXPECT_NEAR(hess[k], htrue[k], 1e-3);
    }
}
