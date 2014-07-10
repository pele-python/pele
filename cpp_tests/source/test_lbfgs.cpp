#include "pele/array.h"
#include "pele/lj.h"
#include "pele/lbfgs.h"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <gtest/gtest.h>
#include <cmath>
#include <memory>

using pele::Array;

TEST(LbfgsLJ, TwoAtom_Works){
    auto lj = std::make_shared<pele::LJ> (1., 1.);
    Array<double> x0(6, 0);
    x0[0] = 2.;
    pele::LBFGS lbfgs(lj, x0);
    lbfgs.run();
    ASSERT_GT(lbfgs.get_nfev(), 1);
    ASSERT_GT(lbfgs.get_niter(), 1);
    ASSERT_LT(lbfgs.get_rms(), 1e-4);
    ASSERT_LT(lbfgs.get_rms(), 1e-4);
    ASSERT_NEAR(lbfgs.get_f(), -.25, 1e-10);
    Array<double> x = lbfgs.get_x();
    double dr, dr2 = 0;
    for (size_t i = 0; i < 3; ++i){
        dr = (x[i] - x[3+i]);
        dr2 += dr * dr;
    }
    dr = sqrt(dr2);
    ASSERT_NEAR(dr, pow(2., 1./6), 1e-5);
    Array<double> g = lbfgs.get_g();
    ASSERT_NEAR(g[0], -g[3], 1e-10);
    ASSERT_NEAR(g[1], -g[4], 1e-10);
    ASSERT_NEAR(g[2], -g[5], 1e-10);
    double rms = pele::norm(g) / sqrt(g.size());
    ASSERT_NEAR(rms, lbfgs.get_rms(), 1e-10);
}

TEST(LbfgsLJ, Reset_Works){
    auto lj = std::make_shared<pele::LJ> (1., 1.);
    Array<double> x0(6, 0);
    x0[0] = 2.;
    // lbfgs1 will minimize straight from x0
    pele::LBFGS lbfgs1(lj, x0);
    lbfgs1.run();

    // lbfgs2 will first minimize from x2 (!=x0) then reset from x0
    // it should end up exactly the same as lbfgs1
    Array<double> x2 = x0.copy();
    x2[1] = 2;
    pele::LBFGS lbfgs2(lj, x2);
    double H0 = lbfgs2.get_H0();
    lbfgs2.run();
    // now reset from x0
    lbfgs2.reset(x0);
    lbfgs2.set_H0(H0);
    lbfgs2.run();

    cout << lbfgs1.get_x() << "\n";
    cout << lbfgs2.get_x() << "\n";

    ASSERT_EQ(lbfgs1.get_nfev(), lbfgs2.get_nfev());
    ASSERT_EQ(lbfgs1.get_niter(), lbfgs2.get_niter());

    for (size_t i=0; i<x0.size(); ++i){
        ASSERT_DOUBLE_EQ(lbfgs1.get_x()[i], lbfgs2.get_x()[i]);
    }
    ASSERT_DOUBLE_EQ(lbfgs1.get_f(), lbfgs2.get_f());
    ASSERT_DOUBLE_EQ(lbfgs1.get_rms(), lbfgs2.get_rms());
//    ASSERT_EQ(lbfgs1.get_niter(), lbfgs1.get_niter());
//    ASSERT_GT(lbfgs.get_niter(), 1);
//    ASSERT_LT(lbfgs.get_rms(), 1e-4);
//    ASSERT_LT(lbfgs.get_rms(), 1e-4);
//    ASSERT_NEAR(lbfgs.get_f(), -.25, 1e-10);


}
