#include <gtest/gtest.h>

#include "pele/lj.h"
#include "pele/steepest_descent.h"

TEST(GDLJ, TwoAtom_Works)
{
    std::shared_ptr<pele::BasePotential> lj = std::make_shared<pele::LJ>(1, 1);
    pele::Array<double> x0(6, 0);
    x0[0] = 2;
    pele::SteepestDescent gd(lj, x0, 1e-1, 1e-5);
    gd.run();
    ASSERT_GT(gd.get_nfev(), 1);
    ASSERT_GT(gd.get_niter(), 1);
    ASSERT_LT(gd.get_rms(), 1e-4);
    ASSERT_NEAR(gd.get_f(), -0.25, 1e-10);
}
