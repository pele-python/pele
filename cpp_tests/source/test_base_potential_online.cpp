#include <gtest/gtest.h>

#include "pele/base_potential.h"
#include "pele/steepest_descent.h"
#include "pele/stochastic_gradient_descent.h"
#include "pele/stochastic_diagonal_levenberg_marquardt.h"
#include "pele/xy_model_online.h"

class BasePotentialOnlineTest : public ::testing::Test {
public:
    double etrue;
    pele::Array<double> x, g, gn, g2, gtrue;
    pele::Array<size_t> head, tail;
    virtual void SetUp()
    {
        head = {0, 0, 1};
        tail = {1, 2, 2};
        x = {M_PI, 0, 0, 0};
        g = {42, 42, 42, 42};
        gn = g.copy();
        g2 = {42, 42, 42, 42};
        gtrue = {0, 0, 0, 0};
        etrue = 2;
    }
};

TEST_F(BasePotentialOnlineTest, Basic_Works)
{
    std::shared_ptr<pele::BasePotentialOnline> pot = std::make_shared<pele::XYModelOnline>(x.size(), head, tail);
    const double e = pot->get_energy(x);
    EXPECT_DOUBLE_EQ(e, etrue);
    EXPECT_DOUBLE_EQ(pot->get_energy(x, 0), 2);
    EXPECT_DOUBLE_EQ(pot->get_energy(x, 1), 0);
    EXPECT_DOUBLE_EQ(pot->get_energy(x, 2), 0);
    EXPECT_DOUBLE_EQ(pot->get_energy(x, 3), 0);
    pot->get_energy_ogradient(x, 0, g);
    pot->numerical_ogradient(x, 0, gn);
    for (size_t i = 0; i < g.size(); ++i) {
        EXPECT_NEAR(g[i], gtrue[i], 1e-10);
        EXPECT_NEAR(gn[i], gtrue[i], 1e-10);
    }
}

class BaseWrapE : public pele::BasePotential {
    std::shared_ptr<pele::BasePotentialOnline> pot;
public:
    BaseWrapE(const size_t N, pele::Array<size_t> h, pele::Array<size_t> t)
        : pot(std::make_shared<pele::XYModelOnline>(N, h, t))
    {}
    double get_energy(pele::Array<double> x)
    {
        return pot->get_energy(x);
    }
};

TEST_F(BasePotentialOnlineTest, XY1D8_Works)
{
    pele::Array<size_t> h8 = {0, 0, 1, 2, 3, 4, 5, 6};
    pele::Array<size_t> t8 = {7, 1, 2, 3, 4, 5, 6, 7};
    std::shared_ptr<pele::BasePotential> pot = std::make_shared<pele::XYModelOnline>(8, h8, t8);
    std::shared_ptr<pele::BasePotential> pot_batch = std::make_shared<BaseWrapE>(8, h8, t8);
    x = pele::Array<double>(8);
    // global minimum
    x.assign(0);
    EXPECT_NEAR(pot->get_energy(x), -16, 1e-10);
    EXPECT_NEAR(pot_batch->get_energy(x), -16, 1e-10);
    EXPECT_NEAR(std::static_pointer_cast<pele::BasePotentialOnline>(pot)->get_energy(x, 3), -2, 1e-10);
    // local minimum
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = i * (M_PI / 4);
    }
    EXPECT_NEAR(pot->get_energy(x), -8 * std::sqrt(2), 1e-10);
    // opimisation: deterministic, stochastic
    x = {0, 1, 2, 3, 4, 5, 6, 7};
    pele::SteepestDescent optimizer_deterministic(pot_batch, x.copy(), 0.1, 1e-8, false);
    pele::StochasticGradientDescent optimizer_stochastic(pot, x.copy(), 0.5, 1e-8, 41, false);
    const double epsilon = 1;
    const double mu = 1;
    const double gamma = 1e-1;
    const double tol = 1e-8;
    const size_t seed = 42;
    const bool verbose = true;
    pele::StochasticDiagonalLevenbergMarquardt sdlm(pot, x.copy(), epsilon, mu, gamma, tol, seed, verbose);
    optimizer_stochastic.set_max_iter(10000);
    optimizer_deterministic.run();
    optimizer_stochastic.run();
    sdlm.run();
    const pele::Array<double> x_deterministic_minimum = optimizer_deterministic.get_x();
    const pele::Array<double> x_stochastic_minimum = optimizer_stochastic.get_x();
    const pele::Array<double> x_sdlm = sdlm.get_x();
    EXPECT_NEAR(pot->get_energy(x_deterministic_minimum), -8 * std::sqrt(2), 1e-10);
    EXPECT_NEAR(pot->get_energy(x_stochastic_minimum), -16, 1e-10);
    EXPECT_NEAR(pot->get_energy(x_sdlm), -16, 1e-10);
    for (size_t i = 1; i < x_deterministic_minimum.size(); ++i) {
        const double delta_deterministic = std::fabs(fmod(x_deterministic_minimum[i], 2 * M_PI) - fmod(x_deterministic_minimum[i - 1], 2 * M_PI));
        EXPECT_NEAR(delta_deterministic, M_PI / 4, 1e-6);
        const double delta_stochastic = std::fabs(fmod(x_stochastic_minimum[i], 2 * M_PI) - fmod(x_stochastic_minimum[i - 1], 2 * M_PI));
        EXPECT_NEAR(delta_stochastic, 0, 1e-6);
    }
}
