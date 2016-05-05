#include <gtest/gtest.h>

#include "pele/xy_model_online.h"

class BasePotentialOnlineTest : public ::testing::Test {
public:
    double etrue;
    pele::Array<double> x, g, g2, gtrue, g2true;
    virtual void SetUp()
    {
        x = {1, 1, 1, 1};
        g = {1, 1, 1, 1};
        g2 = {1, 1, 1, 1};
        gtrue = {1, 1, 1, 1};
        g2true = {0, 0, 0, 0};
        etrue = 42;
    }
};

TEST_F(BasePotentialOnlineTest, OnlyEnergy_Works)
{
    pele::XYModelOnline pot(42);
    /*
    const double e = pot.get_energy(x);
    EXPECT_EQ(1u, pot.call_count);
    EXPECT_NEAR(e, etrue, 1e-10);
    */
}
