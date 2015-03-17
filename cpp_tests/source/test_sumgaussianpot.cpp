#include <gtest/gtest.h>
#include <memory>

#include "pele/array.h"
#include "pele/modified_fire.h"
#include "pele/sumgaussianpot.h"

class TestSumGaussianPot : public ::testing::Test {
public:
    size_t ndim, npot, ndof;
    virtual void SetUp()
    {
        ndim = 2;
        npot = 2;
        ndof = ndim*npot;
    }
};

TEST_F(TestSumGaussianPot, OneGaussWorks) {
    pele::Array<double> mean = pele::Array<double>(ndof, 0);
    pele::Array<double> cov = pele::Array<double>(ndof, 1);
    pele::Array<double> initial_coords = pele::Array<double>(ndof, 1);
    std::shared_ptr<pele::GaussianPot> gauss = std::make_shared<pele::GaussianPot>(mean, cov);
    std::shared_ptr<pele::MODIFIED_FIRE> opt = std::make_shared<pele::MODIFIED_FIRE>(gauss, initial_coords, .1, 1, 1);
    opt->run(1000);
    pele::Array<double> result = opt->get_x();
    for (size_t i = 0; i < ndof; ++i) {
        EXPECT_NEAR(result[i], 0, 1e-4);
    }
    EXPECT_NEAR(gauss->get_energy(mean), -1, 1e-6);
    EXPECT_NEAR(gauss->get_energy(result), -1, 1e-6);
    pele::Array<double> grad(ndof, 42);
    gauss->get_energy_gradient(mean, grad);
    for (size_t i = 0; i < ndof; ++i) {
        EXPECT_DOUBLE_EQ(grad[i], 0.0);
    }
}

TEST_F(TestSumGaussianPot, SumGaussWorks) {
    ndof = 4;
    pele::Array<double> mean = pele::Array<double>(ndof, 0);
    pele::Array<double> cov = pele::Array<double>(ndof, 1);
    pele::Array<double> initial_coords = pele::Array<double>(ndof, 1);
    std::shared_ptr<pele::SumGaussianPot> sumgauss = std::make_shared<pele::SumGaussianPot>(ndim, mean, cov);
    std::shared_ptr<pele::MODIFIED_FIRE> opt = std::make_shared<pele::MODIFIED_FIRE>(sumgauss, initial_coords, .1, 1, 1);
    opt->run(1000);
    pele::Array<double> result = opt->get_x();
    for (size_t i = 0; i < ndim; ++i) {
        EXPECT_NEAR(result[i], 0, 1e-4);
    }
    EXPECT_NEAR(sumgauss->get_energy(mean), -2, 1e-6);
    EXPECT_NEAR(sumgauss->get_energy(result), -2, 1e-6);
    pele::Array<double> grad(ndim, 42);
    sumgauss->get_energy_gradient(mean, grad);
    for (size_t i = 0; i < ndim; ++i) {
        EXPECT_DOUBLE_EQ(grad[i], 0.0);
    }
}

TEST_F(TestSumGaussianPot, SumGaussWorks2) {
    ndof = 4;
    pele::Array<double> mean = pele::Array<double>(ndof, 10);
    pele::Array<double> cov = pele::Array<double>(ndof, 1);
    pele::Array<double> initial_coords = pele::Array<double>(ndof, 9);
    std::shared_ptr<pele::SumGaussianPot> sumgauss = std::make_shared<pele::SumGaussianPot>(ndim, mean, cov);
    std::shared_ptr<pele::MODIFIED_FIRE> opt = std::make_shared<pele::MODIFIED_FIRE>(sumgauss, initial_coords, .1, 1, 1);
    opt->run(1000);
    pele::Array<double> result = opt->get_x();
    for (size_t i = 0; i < ndim; ++i) {
        EXPECT_NEAR(result[i], 10, 1e-4);
    }
    EXPECT_NEAR(sumgauss->get_energy(mean), -2, 1e-6);
    EXPECT_NEAR(sumgauss->get_energy(result), -2, 1e-6);
    pele::Array<double> grad(ndim, 42);
    sumgauss->get_energy_gradient(mean, grad);
    for (size_t i = 0; i < ndim; ++i) {
        EXPECT_DOUBLE_EQ(grad[i], 0.0);
    }
}
