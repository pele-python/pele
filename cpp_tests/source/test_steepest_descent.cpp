#include <cmath>
#include <gtest/gtest.h>

#include "pele/lj.h"
#include "pele/steepest_descent.h"

// see test_lbfgs.cpp

TEST(GDLJ, TwoAtom_Works)
{
    std::shared_ptr<pele::BasePotential> lj = std::make_shared<pele::LJ>(1, 1);
    pele::Array<double> x0(6, 0);
    x0[0] = 2;
    pele::SteepestDescent gd(lj, x0);
    gd.run();
    ASSERT_GT(gd.get_nfev(), 1);
    ASSERT_GT(gd.get_niter(), 1);
    ASSERT_LT(gd.get_rms(), 1e-4);
    ASSERT_NEAR(gd.get_f(), -0.25, 1e-10);
    const pele::Array<double> x = gd.get_x();
    double dr;
    double dr2 = 0;
    for (size_t i = 0; i < 3; ++i) {
        dr = x[i] - x[3 + i];
        dr2 += dr * dr;
    }
    dr = std::sqrt(dr2);
    ASSERT_NEAR(dr, std::pow(2, 1. / 6), 1e-5);
    const pele::Array<double> g = gd.get_g();
    for (size_t i = 0; i < 3; ++i) {
        ASSERT_NEAR(g[i], -g[3 + i], 1e-10);
    }
    const double rms = pele::norm(g) / std::sqrt(g.size());
    ASSERT_NEAR(rms, gd.get_rms(), 1e-10);
}

TEST(GDLJ, Reset_Works)
{
    auto lj = std::make_shared<pele::LJ> (1., 1.);
    pele::Array<double> x0(6, 0);
    x0[0] = 2.;
    // gd1 will minimize straight from x0
    pele::SteepestDescent gd1(lj, x0);
    gd1.run();
    // gd2 will first minimize from x2 (!=x0) then reset from x0
    // it should end up exactly the same as gd1
    pele::Array<double> x2 = x0.copy();
    x2[1] = 2;
    pele::SteepestDescent gd2(lj, x2);
    gd2.run();
    // now reset from x0
    gd2.reset(x0);
    gd2.set_eta(0.1);
    gd2.run();
    ASSERT_EQ(gd1.get_nfev(), gd2.get_nfev());
    ASSERT_EQ(gd1.get_niter(), gd2.get_niter());
    for (size_t i = 0; i < x0.size(); ++i) {
        ASSERT_DOUBLE_EQ(gd1.get_x()[i], gd2.get_x()[i]);
    }
    ASSERT_DOUBLE_EQ(gd1.get_f(), gd2.get_f());
    ASSERT_DOUBLE_EQ(gd1.get_rms(), gd2.get_rms());
}

TEST(GDLJ, SetFuncGradientWorks)
{
    auto lj = std::make_shared<pele::LJ> (1., 1.);
    pele::Array<double> x0(6, 0);
    x0[0] = 2.;
    pele::SteepestDescent gd1(lj, x0);
    pele::SteepestDescent gd2(lj, x0);
    auto grad = x0.copy();
    double e = lj->get_energy_gradient(x0, grad);
    // set the gradient for  gd2.  It should have the same result, but
    // one fewer function evaluation.
    gd2.set_func_gradient(e, grad);
    gd1.run();
    gd2.run();
    ASSERT_EQ(gd1.get_nfev(), gd2.get_nfev() + 1);
    ASSERT_EQ(gd1.get_niter(), gd2.get_niter());
    ASSERT_DOUBLE_EQ(gd1.get_f(), gd2.get_f());
}
