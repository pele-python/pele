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
