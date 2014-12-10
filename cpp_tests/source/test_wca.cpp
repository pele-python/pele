#include "pele/array.h"
#include "pele/wca.h"
#include "pele/hs_wca.h"
#include "pele/meta_pow.h"
#include "test_utils.hpp"

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <gtest/gtest.h>

static double const EPS = std::numeric_limits<double>::min();
#define EXPECT_NEAR_RELATIVE(A, B, T)  EXPECT_NEAR(A/(fabs(A)+fabs(B) + EPS), B/(fabs(A)+fabs(B) + EPS), T)

using pele::Array;
using pele::WCA;
using pele::HS_WCA;
using pele::pos_int_pow;

class WCATest :  public PotentialTest
{
public:
    double sig, eps;
    size_t natoms;

    virtual void setup_potential(){
        pot = std::shared_ptr<pele::BasePotential> (new pele::WCA(sig, eps));
    }

    virtual void SetUp(){
        natoms = 3;
        sig = 1.4;
        eps = 2.1;
        x = Array<double>(natoms*3);
        x[0] = 0.1;
        x[1] = 0.2;
        x[2] = 0.3;
        x[3] = 0.44;
        x[4] = 0.55;
        x[5] = 1.66;
        x[6] = 0.88;
        x[7] = 1.1;
        x[8] = 3.32;
        etrue = 0.9009099166892105;

        setup_potential();
    }
};

TEST_F(WCATest, Energy_Works){
    test_energy();
}
TEST_F(WCATest, EnergyGradient_AgreesWithNumerical){
    test_energy_gradient();
}
TEST_F(WCATest, EnergyGradientHessian_AgreesWithNumerical){
    test_energy_gradient_hessian();
}

class WCAAtomListTest :  public WCATest
{
public:
    virtual void setup_potential(){
        pele::Array<size_t> atoms(natoms);
        for (size_t i =0; i<atoms.size(); ++i){
            atoms[i] = i;
        }
        pot = std::shared_ptr<pele::BasePotential> (new pele::WCAAtomList(
                sig, eps, atoms
                ));
    }
};

TEST_F(WCAAtomListTest, Energy_Works){
    test_energy();
}
TEST_F(WCAAtomListTest, EnergyGradient_AgreesWithNumerical){
    test_energy_gradient();
}
TEST_F(WCAAtomListTest, EnergyGradientHessian_AgreesWithNumerical){
    test_energy_gradient_hessian();
}



/*
 * HS_WCA tests
 */
class HS_WCATest :  public ::testing::Test
{
public:
    double eps, sca, etrue;
    Array<double> x, g, gnum, radii;
    virtual void SetUp(){
        eps = 2.1;
        sca = 1.4;
        x = Array<double>(9);
        x[0] = 0.1;
        x[1] = 0.2;
        x[2] = 0.3;
        x[3] = 0.44;
        x[4] = 0.55;
        x[5] = 1.66;
        x[6] = 0.88;
        x[7] = 1.1;
        x[8] = 3.32;
        radii = Array<double>(3);
        double f = .35;
        radii[0] = .91 * f;
        radii[1] = 1.1 * f;
        radii[2] = 1.13 * f;
        etrue = 189.41835811474974;
        g = Array<double>(x.size());
        gnum = Array<double>(x.size());
    }
};

TEST_F(HS_WCATest, Energy_Works){
    HS_WCA<3> pot(eps, sca, radii);
    double e = pot.get_energy(x);
    ASSERT_NEAR(e, etrue, 1e-10);
}

class OtherfHS_WCA {
public:
    OtherfHS_WCA(const double r_sum_, const double infinity_, const double epsilon_, const double alpha_)
        : r_sum(r_sum_),
          infinity(infinity_),
          epsilon(epsilon_),
          alpha(alpha_),
          r_sum_soft((1 + alpha) * r_sum)
    {}
    double operator()(const double r) const
    {
        if (r >= r_sum_soft) {
            return 0;
        }
        if (r <= r_sum) {
            return infinity;
        }
        const double numerator = sigma();
        const double denominator = pos_int_pow<2>(r) - pos_int_pow<2>(r_sum);
        const double ratio = numerator / denominator;
        return std::max<double>(0, epsilon * (4 * (pos_int_pow<12>(ratio) - pos_int_pow<6>(ratio)) + 1));
    }
    double grad(const double r) const
    {
        if (r >= r_sum_soft) {
            return 0;
        }
        if (r <= r_sum) {
            return 0;
        }
        return (-8 * epsilon * r * sigma()) / pos_int_pow<2>(pos_int_pow<2>(r) - pos_int_pow<2>(r_sum)) * (12 * pos_int_pow<11>(g(r)) - 6 * pos_int_pow<5>(g(r)));
    }
    // This is not the gradient. Rather grad(r) / (-r).
    double scaled_grad(const double r) const
    {
        return grad(r) / (-r);
    }
    double sigma() const { return (2 * alpha + pos_int_pow<2>(alpha)) * pos_int_pow<2>(r_sum) * std::pow(2, -static_cast<double>(1) / static_cast<double>(6)); }
    double g(const double r) const { return sigma() / (pos_int_pow<2>(r) - pos_int_pow<2>(r_sum)); }
private:
    const double r_sum;
    const double infinity;
    const double epsilon;
    const double alpha;
    const double r_sum_soft;
};

class OthersfHS_WCA {
public:
    OthersfHS_WCA(const double r_sum_, const double epsilon_, const double alpha_, const double delta_=1e-10)
        : r_sum(r_sum_),
          epsilon(epsilon_),
          alpha(alpha_),
          delta(delta_),
          fHS_WCA(r_sum, std::numeric_limits<double>::max(), epsilon, alpha),
          r_sum_soft((1 + alpha) * r_sum),
          r_X(r_sum + delta)
    {}
    double operator()(const double r) const
    {
        if (r > r_sum_soft) {
            return 0;
        }
        if (r > r_X) {
            return fHS_WCA(r);
        }
        //return fHS_WCA(r_X) - (r - r_X) * fHS_WCA.scaled_grad(r_X);
        //return fHS_WCA(r_X) - (r - r_X) * fHS_WCA.grad(r_X) / (-r_X);
        return fHS_WCA(r_X) + (r - r_X) * fHS_WCA.grad(r_X);
    }
private:
    const double r_sum;
    const double epsilon;
    const double alpha;
    const double delta;
    const OtherfHS_WCA fHS_WCA;
    const double r_sum_soft;
    const double r_X;
};

TEST_F(HS_WCATest, ExtendedEnergyTest_Works){
    HS_WCA<3> pot(eps, sca, radii);
    const double e = pot.get_energy(x);
    EXPECT_DOUBLE_EQ(e, etrue);
    pele::HS_WCA_interaction pair_pot(eps, sca, radii);
    const size_t atom_a = 0;
    const size_t atom_b = 1;
    const double r_sum = radii[atom_a] + radii[atom_b];
    const double rmin = r_sum / 2;
    const size_t nr_points = 10000;
    const double rmax = 3 * r_sum * (sca + 1);
    const double rdelta = (rmax - rmin) / (nr_points - 1);
    const double infinity = pair_pot._infty;
    OtherfHS_WCA other_implementation(r_sum, infinity, eps, sca);
    for (size_t i = 0; i < nr_points; ++i) {
        const double r = rmin + i * rdelta;
        const double e_pair_pot_f_energy = pair_pot.energy(pos_int_pow<2>(r), atom_a, atom_b);
        double gab;
        const double e_pair_pot_f_energy_gradient = pair_pot.energy_gradient(pos_int_pow<2>(r), &gab, atom_a, atom_b);
        double hab;
        const double e_pair_pot_f_energy_gradient_hessian = pair_pot.energy_gradient_hessian(pos_int_pow<2>(r), &gab, &hab, atom_a, atom_b);
        EXPECT_LE(0, e_pair_pot_f_energy);
        EXPECT_DOUBLE_EQ(e_pair_pot_f_energy, e_pair_pot_f_energy_gradient);
        EXPECT_LE(0, e_pair_pot_f_energy_gradient);
        EXPECT_DOUBLE_EQ(e_pair_pot_f_energy, e_pair_pot_f_energy_gradient_hessian);
        EXPECT_LE(0, e_pair_pot_f_energy_gradient_hessian);
        const double e_other = other_implementation(r);
        EXPECT_NEAR_RELATIVE(e_pair_pot_f_energy, e_other, 1e-10);
        EXPECT_LE(0, e_other);
        if (r > (sca + 1) * r_sum) {
            EXPECT_DOUBLE_EQ(e_other, 0);
            EXPECT_DOUBLE_EQ(e_pair_pot_f_energy, 0);
        }
        else if (r < r_sum) {
            EXPECT_DOUBLE_EQ(e_other, infinity);
            EXPECT_DOUBLE_EQ(e_pair_pot_f_energy, infinity);
        }
    }
    // Un-comment the following to plot sfHS-WCA potential.
    /*
    std::vector<double> x_;
    std::vector<double> y_;
    const size_t scale = 2;
    //pele::sf_HS_WCA_interaction(eps, sca, radii).evaluate_pair_potential(rmin, rmax, nr_points, atom_a, atom_b, x_, y_);
    pele::sf_HS_WCA_interaction(eps, sca, radii).evaluate_pair_potential(0.70348, 0.70352, scale * nr_points, atom_a, atom_b, x_, y_);
    //std::ofstream out("test_sfhs_wca_shape.txt");
    std::ofstream out("test_sfhs_wca_shape_zoom_grad.txt");
    out.precision(std::numeric_limits<double>::digits10);
    for (size_t i = 0; i < scale * nr_points; ++i) {
        out << x_.at(i) << "\t" << y_.at(i) << "\n";
    }
    out.close();
    */
    // Below HS_WCA_interaction's infinity, HS_WCA_interaction and
    // sf_HS_WCA_interaction have to be the same.
    // Also sf_HS_WCA_interaction should agree with the second
    // alternative implementation given above, for all points.
    pele::sf_HS_WCA_interaction sf_pair_pot(eps, sca, radii);
    OthersfHS_WCA sf_other_implementation(r_sum, eps, sca);
    for (size_t i = 0; i < nr_points; ++i) {
        const double r = rmin + i * rdelta;
        double pair_pot_gab;
        const double pair_pot_e = pair_pot.energy_gradient(pos_int_pow<2>(r), &pair_pot_gab, atom_a, atom_b);
        double sf_pair_pot_gab;
        const double sf_pair_pot_e = sf_pair_pot.energy_gradient(pos_int_pow<2>(r), &sf_pair_pot_gab, atom_a, atom_b);
        if (pair_pot_e < infinity) {
            // Here sf and f have to give the same result.
            EXPECT_LE(0, pair_pot_e);
            EXPECT_LE(0, sf_pair_pot_e);
            EXPECT_DOUBLE_EQ(pair_pot_e, sf_pair_pot_e);
            EXPECT_DOUBLE_EQ(pair_pot_gab, sf_pair_pot_gab);
        }
        // sf always has to agree with a simpler implementation given above.
        // Numerical differences can be present.
        const double alternative_sf_pair_pot_e = sf_other_implementation(r);
        EXPECT_LE(0, alternative_sf_pair_pot_e);
        EXPECT_NEAR_RELATIVE(sf_pair_pot_e, alternative_sf_pair_pot_e, 1e-10);
    }
}

TEST_F(HS_WCATest, EnergyGradient_AgreesWithNumerical){
    HS_WCA<3> pot(eps, sca, radii);
    double e = pot.get_energy_gradient(x, g);
    double ecomp = pot.get_energy(x);
    ASSERT_NEAR(e, ecomp, 1e-10);
    pot.numerical_gradient(x, gnum, 1e-6);
    for (size_t k=0; k<6; ++k){
        ASSERT_NEAR(g[k], gnum[k], 1e-6);
    }
}

TEST_F(HS_WCATest, EnergyGradientHessian_AgreesWithNumerical){
    HS_WCA<3> pot(eps, sca, radii);
    Array<double> h(x.size()*x.size());
    Array<double> hnum(h.size());
    double e = pot.get_energy_gradient_hessian(x, g, h);
    double ecomp = pot.get_energy(x);
    pot.numerical_gradient(x, gnum);
    pot.numerical_hessian(x, hnum);

    EXPECT_NEAR(e, ecomp, 1e-10);
    for (size_t i=0; i<g.size(); ++i){
        ASSERT_NEAR(g[i], gnum[i], 1e-6);
    }
    for (size_t i=0; i<h.size(); ++i){
        ASSERT_NEAR(h[i], hnum[i], 1e-3);
    }
}
