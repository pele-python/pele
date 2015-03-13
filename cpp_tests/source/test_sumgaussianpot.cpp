#include <gtest/gtest.h>

#include "pele/array.h"
#include "pele/modified_fire.h"
#include "pele/sumgaussianpot.h"

class TestSumGaussianPot : public ::testing::Test {
public:
    size_t ndim;
    std::shared_ptr<pele::SumGaussianPot> pot;
    std::shared_ptr<pele::MODIFIED_FIRE> opt;
    pele::Array<double> correct;
    pele::Array<double> ini;
    virtual void SetUp()
    {
        ndim = 40;
        correct = pele::Array<double>(ndim, 1);
        ini = pele::Array<double>(ndim, 42);
        pot = std::make_shared<pele::SumGaussianPot>(1, correct, correct);
        opt = std::make_shared<pele::MODIFIED_FIRE>(pot, ini, .1, 1, 1);
    }
};

TEST_F(TestSumGaussianPot, Works) {
    opt->reset(ini);
    opt->run();
    pele::Array<double> result = opt->get_x();
    EXPECT_DOUBLE_EQ(pot->get_energy(correct), -static_cast<double>(ndim));
    /*
    for (size_t i = 0; i < ndim; ++i) {
        EXPECT_NEAR(result[i], correct[i], 1e-5);
    }
    */
}
