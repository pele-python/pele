#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include "pele/lj.h"
#include "pele/modified_fire.h"

using pele::Array;

TEST(FireLJ, TwoAtom_Works){
    auto lj = std::make_shared<pele::LJ> (1., 1.);
    Array<double> x0(6, 0);
    x0[0] = 2.;
    pele::MODIFIED_FIRE fire(lj, x0, 1, 1, 1);
    fire.set_iprint(1);
    fire.run();
    ASSERT_GT(fire.get_nfev(), 1);
    ASSERT_GT(fire.get_niter(), 1);
    ASSERT_LT(fire.get_rms(), 1e-4);
    ASSERT_LT(fire.get_rms(), 1e-4);
    ASSERT_NEAR(fire.get_f(), -.25, 1e-8);
    Array<double> x = fire.get_x();
    double dr, dr2 = 0;
    for (size_t i = 0; i < 3; ++i){
        dr = (x[i] - x[3+i]);
        dr2 += dr * dr;
    }
    dr = sqrt(dr2);
    ASSERT_NEAR(dr, pow(2., 1./6), 1e-5);
    Array<double> g = fire.get_g();
    ASSERT_NEAR(g[0], -g[3], 1e-10);
    ASSERT_NEAR(g[1], -g[4], 1e-10);
    ASSERT_NEAR(g[2], -g[5], 1e-10);
    double rms = pele::norm(g) / sqrt(g.size());
    ASSERT_NEAR(rms, fire.get_rms(), 1e-10);
}

TEST(FireLJ, Reset_Works){
    auto lj = std::make_shared<pele::LJ> (1., 1.);
    Array<double> x0(6, 0);
    x0[0] = 2.;
    // lbfgs1 will minimize straight from x0
    pele::MODIFIED_FIRE fire1(lj, x0, 1, 1, 1);
    fire1.run();

    // lbfgs2 will first minimize from x2 (!=x0) then reset from x0
    // it should end up exactly the same as lbfgs1
    Array<double> x2 = x0.copy();
    x2[1] = 2;
    pele::MODIFIED_FIRE fire2(lj, x2, 1, 1, 1);
    fire2.run();
    // now reset from x0
    fire2.reset(x0);
    fire2.run();

    std::cout << fire1.get_x() << "\n";
    std::cout << fire2.get_x() << "\n";

    ASSERT_EQ(fire1.get_nfev(), fire2.get_nfev());
    ASSERT_EQ(fire1.get_niter(), fire2.get_niter());

    for (size_t i=0; i<x0.size(); ++i){
        ASSERT_DOUBLE_EQ(fire1.get_x()[i], fire2.get_x()[i]);
    }
    ASSERT_DOUBLE_EQ(fire1.get_f(), fire2.get_f());
    ASSERT_DOUBLE_EQ(fire1.get_rms(), fire2.get_rms());
//    ASSERT_EQ(lbfgs1.get_niter(), lbfgs1.get_niter());
//    ASSERT_GT(lbfgs.get_niter(), 1);
//    ASSERT_LT(lbfgs.get_rms(), 1e-4);
//    ASSERT_LT(lbfgs.get_rms(), 1e-4);
//    ASSERT_NEAR(lbfgs.get_f(), -.25, 1e-10);


}


TEST(FireLJ, SetFuncGradientWorks){
    auto lj = std::make_shared<pele::LJ> (1., 1.);
    Array<double> x0(6, 0);
    x0[0] = 2.;
    pele::MODIFIED_FIRE fire1(lj, x0, 1, 1, 1);
    pele::MODIFIED_FIRE fire2(lj, x0, 1, 1, 1);
    auto grad = x0.copy();
    double e = lj->get_energy_gradient(x0, grad);

    // set the gradient for  lbfgs2.  It should have the same result, but
    // one fewer function evaluation.
    fire2.set_func_gradient(e, grad);
    fire1.run();
    fire2.run();
    ASSERT_EQ(fire1.get_nfev(), fire2.get_nfev() + 1);
    ASSERT_EQ(fire1.get_niter(), fire2.get_niter());
    ASSERT_DOUBLE_EQ(fire1.get_f(), fire2.get_f());
}


