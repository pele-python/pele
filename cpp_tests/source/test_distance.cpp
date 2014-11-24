#include <random>
#include <numeric>

#include <gtest/gtest.h>

#include "pele/distance.h"

using pele::Array;
using pele::periodic_distance;
using pele::cartesian_distance;

static double const EPS = std::numeric_limits<double>::min();
#define EXPECT_NEAR_RELATIVE(A, B, T)  ASSERT_NEAR(A/(fabs(A)+fabs(B) + EPS), B/(fabs(A)+fabs(B) + EPS), T)

class DistanceTest :  public ::testing::Test {
public:
    double* x2;
    double* x3;
    double* x42;
    double* y2;
    double* y3;
    double* y42;
    virtual void SetUp()
    {
        std::mt19937_64 gen(42);
        std::uniform_real_distribution<double> dist(1, 10);
        x2 = new double[2];
        y2 = new double[2];
        x3 = new double[3];
        y3 = new double[3];
        x42 = new double[42];
        y42 = new double[42];
        for (size_t i = 0; i < 2; ++i) {
            x2[i] = dist(gen);
            y2[i] = dist(gen);
        }
        for (size_t i = 0; i < 3; ++i) {
            x3[i] = dist(gen);
            y3[i] = dist(gen);
        }
        for (size_t i = 0; i < 42; ++i) {
            x42[i] = dist(gen);
            y42[i] = dist(gen);
        }
    }
    virtual void TearDown()
    {
        delete[] x2;
        delete[] y2;
        delete[] x3;
        delete[] y3;
        delete[] x42;
        delete[] y42;
    }
};


TEST_F(DistanceTest, CartesianDistanceNorm_Works)
{
    // compute with pele
    double dx_p_2[2];
    double dx_p_3[3];
    double dx_p_42[42];
    cartesian_distance<2>().get_rij(dx_p_2, x2, y2);
    cartesian_distance<3>().get_rij(dx_p_3, x3, y3);
    cartesian_distance<42>().get_rij(dx_p_42, x42, y42);
    double ds_p_2 = 0;
    double ds_p_3 = 0;
    double ds_p_42 = 0;
    // compute with std
    double dx2[2];
    double dx3[3];
    double dx42[42];
    for (size_t i = 0; i < 2; ++i) {
        dx2[i] = x2[i] - y2[i];
        ASSERT_DOUBLE_EQ(dx_p_2[i], dx2[i]);
        ds_p_2 += dx_p_2[i] * dx_p_2[i];
    }
    for (size_t i = 0; i < 3; ++i) {
        dx3[i] = x3[i] - y3[i];
        ASSERT_DOUBLE_EQ(dx_p_3[i], dx3[i]);
        ds_p_3 += dx_p_3[i] * dx_p_3[i];
    }
    for (size_t i = 0; i < 42; ++i) {
        dx42[i] = x42[i] - y42[i];
        ASSERT_DOUBLE_EQ(dx_p_42[i], dx42[i]);
        ds_p_42 += dx_p_42[i] * dx_p_42[i];
    }
    const double ds2 = std::inner_product(dx2, dx2 + 2, dx2, double(0));
    const double ds3 = std::inner_product(dx3, dx3 + 3, dx3, double(0));
    const double ds42 = std::inner_product(dx42, dx42 + 42, dx42, double(0));
    // compare norms
    ASSERT_DOUBLE_EQ(ds_p_2, ds2);
    ASSERT_DOUBLE_EQ(ds_p_3, ds3);
    ASSERT_DOUBLE_EQ(ds_p_42, ds42);
}

TEST_F(DistanceTest, SimplePeriodicNorm_Works)
{
    // compute with pele
    double dx_p_2[2];
    double dx_p_3[3];
    double dx_p_42[42];
    pele::Array<double> bv2(2, 18);
    pele::Array<double> bv3(3, 18);
    pele::Array<double> bv42(42, 18);
    periodic_distance<2>(bv2).get_rij(dx_p_2, x2, y2);
    periodic_distance<3>(bv3).get_rij(dx_p_3, x3, y3);
    periodic_distance<42>(bv42).get_rij(dx_p_42, x42, y42);
    double ds_p_2 = 0;
    double ds_p_3 = 0;
    double ds_p_42 = 0;
    // compute with std
    double dx2[2];
    double dx3[3];
    double dx42[42];
    for (size_t i = 0; i < 2; ++i) {
        dx2[i] = x2[i] - y2[i];
        ASSERT_DOUBLE_EQ(dx_p_2[i], dx2[i]);
        ds_p_2 += dx_p_2[i] * dx_p_2[i];
    }
    for (size_t i = 0; i < 3; ++i) {
        dx3[i] = x3[i] - y3[i];
        ASSERT_DOUBLE_EQ(dx_p_3[i], dx3[i]);
        ds_p_3 += dx_p_3[i] * dx_p_3[i];
    }
    for (size_t i = 0; i < 42; ++i) {
        dx42[i] = x42[i] - y42[i];
        ASSERT_DOUBLE_EQ(dx_p_42[i], dx42[i]);
        ds_p_42 += dx_p_42[i] * dx_p_42[i];
    }
    const double ds2 = std::inner_product(dx2, dx2 + 2, dx2, double(0));
    const double ds3 = std::inner_product(dx3, dx3 + 3, dx3, double(0));
    const double ds42 = std::inner_product(dx42, dx42 + 42, dx42, double(0));
    // compare norms
    ASSERT_DOUBLE_EQ(ds_p_2, ds2);
    ASSERT_DOUBLE_EQ(ds_p_3, ds3);
    ASSERT_DOUBLE_EQ(ds_p_42, ds42);
}

TEST_F(DistanceTest, NearestImageConvention_Works)
{
    double x_out_of_box2[2];
    double x_boxed_true2[2];
    double x_boxed_per2[2];
    double x_out_of_box3[3];
    double x_boxed_true3[3];
    double x_boxed_per3[3];
    double x_out_of_box42[42];
    double x_boxed_per42[42];
    // The following means that "in box" is in [-8, 8].
    const double L = 16;
    x_out_of_box2[0] = -10;
    x_out_of_box2[1] = 20;
    x_boxed_true2[0] = 6;
    x_boxed_true2[1] = 4;
    std::copy(x_out_of_box2, x_out_of_box2 + 2, x_boxed_per2);
    pele::Array<double> xp2(x_boxed_per2, 2);
    periodic_distance<2>(pele::Array<double>(2, L)).put_in_box(xp2);
    for (size_t i = 0; i < 2; ++i) {
        EXPECT_DOUBLE_EQ(x_boxed_per2[i], x_boxed_true2[i]);
    }
    x_out_of_box3[0] = -9;
    x_out_of_box3[1] = 8.25;
    x_out_of_box3[2] = 12.12;
    x_boxed_true3[0] = 7;
    x_boxed_true3[1] = -7.75;
    x_boxed_true3[2] = -3.88;
    std::copy(x_out_of_box3, x_out_of_box3 + 3, x_boxed_per3);
    pele::Array<double> xp3(x_boxed_per3, 3);
    periodic_distance<3>(pele::Array<double>(3, L)).put_in_box(xp3);
    for (size_t i = 0; i < 3; ++i) {
        EXPECT_DOUBLE_EQ(x_boxed_per3[i], x_boxed_true3[i]);
    }
    // Assert that putting in box is irrelevant for distances.
    std::mt19937_64 gen(42);
    std::uniform_real_distribution<double> dist(-100, 100);
    for (size_t i = 0; i < 42; ++i) {
        x_out_of_box42[i] = dist(gen);
    }
    std::vector<double> ones(42, 1);
    double delta42[42];
    periodic_distance<42>(pele::Array<double>(42, L)).get_rij(delta42, &*ones.begin(), x_out_of_box42);
    const double d2_42_before = std::inner_product(delta42, delta42 + 42, delta42, double(0));
    std::copy(x_out_of_box42, x_out_of_box42 + 42, x_boxed_per42);
    pele::Array<double> xp42(x_boxed_per42, 42);
    periodic_distance<42>(pele::Array<double>(42, L)).put_in_box(xp42);
    for (size_t i = 0; i < 42; ++i) {
        EXPECT_LE(x_boxed_per42[i], 0.5 * L);
        EXPECT_LE(-0.5 * L, x_boxed_per42[i]);
    }
    periodic_distance<42>(pele::Array<double>(42, L)).get_rij(delta42, &*ones.begin(), x_boxed_per42);
    const double d2_42_after = std::inner_product(delta42, delta42 + 42, delta42, double(0));
    EXPECT_DOUBLE_EQ(d2_42_before, d2_42_after);
    // Test for multiple particles.
    const double d2_42_ndim2_before = d2_42_before;
    periodic_distance<2>(pele::Array<double>(2, L)).get_rij(delta42, &*ones.begin(), x_boxed_per42);
    const double d2_42_ndim2_after = std::inner_product(delta42, delta42 + 42, delta42, double(0));
    EXPECT_DOUBLE_EQ(d2_42_ndim2_before, d2_42_ndim2_after);
    std::copy(x_out_of_box42, x_out_of_box42 + 42, x_boxed_per42);
    periodic_distance<2>(pele::Array<double>(2, L)).put_in_box(xp42);
    for (size_t i = 0; i < 42; ++i) {
        EXPECT_LE(x_boxed_per42[i], 0.5 * L);
        EXPECT_LE(-0.5 * L, x_boxed_per42[i]);
    }
}


