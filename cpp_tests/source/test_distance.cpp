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


