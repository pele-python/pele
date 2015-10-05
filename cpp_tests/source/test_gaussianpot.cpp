#include <gtest/gtest.h>
#include <memory>

#include "pele/array.h"
#include "pele/modified_fire.h"
#include "pele/gaussianpot.h"

class TestGaussianPot : public ::testing::Test {
public:
    size_t ndim, npot, ndof;
    virtual void SetUp()
    {
        ndim = 2;
        npot = 2;
        ndof = ndim * npot;
    }
};

TEST_F(TestGaussianPot, OneGaussWorks)
{
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

void energy_test(std::initializer_list<double> x_, const double e_true)
{
    pele::Array<double> x(x_);
    pele::Array<double> mean(x.size(), 0);
    pele::Array<double> cov(x.size(), 1);
    pele::GaussianPot gauss(mean, cov);
    EXPECT_DOUBLE_EQ(e_true, gauss.get_energy(x));
}

void gradient_test(std::initializer_list<double> x_, std::initializer_list<double> grad_true_)
{
    pele::Array<double> x(x_);
    pele::Array<double> grad_true(grad_true_);
    pele::Array<double> mean(x.size(), 0);
    pele::Array<double> cov(x.size(), 1);
    pele::Array<double> grad(x.size());
    pele::GaussianPot(mean, cov).get_energy_gradient(x, grad);
    for (unsigned int i = 0; i < grad.size(); ++i) {
        EXPECT_DOUBLE_EQ(grad_true[i], grad[i]);
    }
}

TEST_F(TestGaussianPot, OneGaussWorksNonZero)
{
    energy_test({0, 0, 0, 0}, -1);
    energy_test({1, 1, 1, 1}, -0.1353352832366127);
    energy_test({2, 2, 2, 2}, -0.0003354626279025118);
    energy_test({1, 2, 3, 4}, -3.059023205018258e-7);
    gradient_test({0, 0, 0, 0}, {0, 0, 0, 0});
    gradient_test({0, 0.1, 0.2, 0.3}, {0, 0.09323938199059482, 0.1864787639811896, 0.2797181459717845});
    gradient_test({1, 1, 1, 1}, {0.1353352832366127, 0.1353352832366127, 0.1353352832366127, 0.1353352832366127});
}
