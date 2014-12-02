#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "pele/harmonic.h"
#include "pele/meta_pow.h"

#include "test_utils.hpp"

using pele::Array;
using pele::pos_int_pow;

static double const EPS = std::numeric_limits<double>::min();
#define EXPECT_NEAR_RELATIVE(A, B, T)  EXPECT_NEAR(A/(fabs(A)+fabs(B) + EPS), B/(fabs(A)+fabs(B) + EPS), T)

class HarmonicTest : public ::testing::Test {
public:
    size_t nr_particles;
    size_t box_dimension;
    size_t nr_dof;
    Array<double> x;
    double k;
    double displ;
    Array<double> xd;
    double e_ini_true;
    void SetUp()
    {
        nr_particles = 100;
        box_dimension = 3;
        nr_dof = nr_particles * box_dimension;
        x = Array<double>(nr_dof, 0);
        k = 12.12;
        displ = 43.44;
        xd = Array<double>(nr_dof, displ);
        e_ini_true = 0.5 * k * nr_dof * pos_int_pow<2>(displ);
    }
};

TEST_F(HarmonicTest, Works) {
    auto pot = std::make_shared<pele::Harmonic>(x, k, box_dimension);
    const double e_ini = pot->get_energy(xd);
    EXPECT_NEAR_RELATIVE(e_ini, e_ini_true, 1e-10);
    for (size_t k = 0; k < 100; ++k) {
        pot->set_k(k);
        const double e_pot = pot->get_energy(xd);
        const double e_true = 0.5 * k * nr_dof * pos_int_pow<2>(displ);
        EXPECT_NEAR_RELATIVE(e_pot, e_true, 1e-10);
        Array<double> actual_grad(nr_dof, 0);
        const double e_pot_grad = pot->get_energy_gradient(xd, actual_grad);
        EXPECT_NEAR_RELATIVE(e_pot, e_pot_grad, 1e-10);
        for (size_t i = 0; i < nr_dof; ++i) {
            EXPECT_NEAR_RELATIVE(actual_grad[i], k * displ, 1e-10);
        }
    }
}

class HarmonicAtomListTest :  public PotentialTest {
public:
    double natoms;
    double k;

    virtual void setup_potential(){
        pele::Array<size_t> atoms(natoms);
        for (size_t i =0; i<atoms.size(); ++i){
            atoms[i] = i;
        }
        pot = std::shared_ptr<pele::BasePotential> (new pele::HarmonicAtomList(
                k, atoms
                ));
    }


    virtual void SetUp()
    {
        natoms = 3;
        k = 1.;
        x = Array<double>(3*natoms);
        x[0]  = 0.1;
        x[1]  = 0.2;
        x[2]  = 0.3;
        x[3]  = 0.44;
        x[4]  = 0.55;
        x[5]  = 1.66;

        x[6] = 0;
        x[7] = 0;
        x[8] = -3.;

        etrue = 17.6197;

        setup_potential();
    }
};

TEST_F(HarmonicAtomListTest, Energy_Works){
    test_energy();
}

TEST_F(HarmonicAtomListTest, EnergyGradient_AgreesWithNumerical){
    test_energy_gradient();
}

TEST_F(HarmonicAtomListTest, EnergyGradientHessian_AgreesWithNumerical){
    test_energy_gradient_hessian();
}

class HarmonicNeighborListTest :  public HarmonicAtomListTest {
public:
    virtual void setup_potential(){
        pele::Array<size_t> nlist(natoms * (natoms - 1) );
        size_t count = 0;
        for (size_t i = 0; i < natoms; ++i){
            for (size_t j = i+1; j < natoms; ++j){
                nlist[count++] = i;
                nlist[count++] = j;
            }
        }
        assert(count == nlist.size());
        pot = std::shared_ptr<pele::BasePotential> (new pele::HarmonicNeighborList(
                k, nlist));
    }
};

TEST_F(HarmonicNeighborListTest, Energy_Works){
    test_energy();
}

TEST_F(HarmonicNeighborListTest, EnergyGradient_AgreesWithNumerical){
    test_energy_gradient();
}

TEST_F(HarmonicNeighborListTest, EnergyGradientHessian_AgreesWithNumerical){
    test_energy_gradient_hessian();
}

