#include <cmath>

#include <gtest/gtest.h>

#include "pele/inversepower.h"
#include "pele/pressure_tensor.h"

TEST(Pressure, BasicComponents_Work)
{
    const size_t nr_dim = 3;
    const size_t power = 2;
    const double eps = 1;
    const size_t nr_particles = 2;
    const double radius = 1;
    pele::Array<double> radii(nr_particles, radius);
    std::shared_ptr<pele::BasePotential> pot = std::make_shared<pele::InverseIntPower<nr_dim, power> >(eps, radii);
    pele::Array<double> x(nr_dim * nr_particles);
    for (size_t i = 0; i < nr_dim; ++i) {
        x[i] = 1;
        x[i + nr_dim] = 1.5;
    }
    const double e = pot->get_energy(x);
    pele::Array<double> delta_x(nr_dim);
    for (size_t i = 0; i < nr_dim; ++i) {
        delta_x[i] = x[i] - x[i + nr_dim];
    }
    const double r = pele::norm(delta_x);
    const double e_true = r < 2 * radius ? eps / power * std::pow((1 - r / (2 * radius)), power) : 0;
    EXPECT_DOUBLE_EQ(e_true, e);
    double dr[nr_dim];
    dynamic_cast<pele::SimplePairwisePotentialInterface*>(pot.get())->get_rij(dr, x.data(), x.data() + nr_dim);
    for (size_t i = 0; i < nr_dim; ++i) {
        EXPECT_DOUBLE_EQ(delta_x[i], dr[i]);
    }
    double gij;
    const double e_i = dynamic_cast<pele::SimplePairwisePotentialInterface*>(pot.get())->get_interaction_energy_gradient(pele::pos_int_pow<2>(r), &gij, 0, 1);
    EXPECT_DOUBLE_EQ(e_true, e_i);
    const double sigma = 2 * radius;
    const double gij_true = eps / (sigma * r) * std::pow(1 - r / sigma, power - 1);
    pele::Array<double> ptensor(nr_dim * nr_dim);
    EXPECT_DOUBLE_EQ(gij_true, gij);
    const double volume = 42;
    const double p = pele::pressure_tensor(pot, x, ptensor, volume);
    pele::Array<double> half_delta_x = delta_x;
    pele::Array<double> force = gij_true * delta_x;
    pele::Array<double> ptensor_true(nr_dim * nr_dim, 0);
    for (size_t k = 0; k < nr_dim; ++k) {
        for (size_t l = k; l < nr_dim; ++l) {
            ptensor_true[k * nr_dim + l] = half_delta_x[k] * force[l];
            ptensor_true[l * nr_dim + k] = half_delta_x[l] * force[k];
        }
    }
    ptensor_true /= volume;
    double tr = 0;
    for (size_t i = 0; i < nr_dim; ++i) {
        tr += ptensor_true[i * nr_dim + i];
    }
    const double p_true = tr / nr_dim;
    EXPECT_DOUBLE_EQ(p_true, p);
    for (size_t i = 0; i < ptensor.size(); ++i) {
        //std::cout << "ptensor[i] / ptensor_true[i]: " << ptensor[i] / ptensor_true[i] << "\n";
        EXPECT_DOUBLE_EQ(ptensor[i], ptensor_true[i]);
    }
}

TEST(Pressure, Interface_Throws)
{
    pele::SimplePairwisePotentialInterface i;
    EXPECT_THROW(i.get_ndim(), std::runtime_error);
    EXPECT_THROW(i.get_rij(NULL, NULL, NULL), std::runtime_error);
    EXPECT_THROW(i.get_interaction_energy_gradient(42.42, NULL, 42, 42), std::runtime_error);
}
