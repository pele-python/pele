#include "pele/array.h"
#include "pele/lj.h"
#include "pele/lj_cut.h"
#include "test_utils.hpp"

#include <iostream>
#include <stdexcept>
#include <vector>
#include <gtest/gtest.h>
#include <cmath>
#include <memory>

using pele::Array;


TEST(LJInteractionTest, Energy_Works){
    double r2 = 1.1;
    double c6 = 1.2;
    double c12 = 2.3;
    pele::lj_interaction ljint(c6, c12);
    ASSERT_NEAR(ljint.energy(r2, 1, 2), 0.39671227804179443, 1e-10);
}

TEST(LJInteractionTest, EnergyGradient_Works){
    double r2 = 1.1;
    double c6 = 1.2;
    double c12 = 2.3;
    pele::lj_interaction ljint(c6, c12);
    double g = 0;
    double e = ljint.energy_gradient(r2, &g, 1, 2);
    ASSERT_NEAR(e, 0.39671227804179443, 1e-10);
    ASSERT_NEAR(g, 9.2454671845389917, 1e-10);
}

TEST(LJInteractionTest, Hessian_Works){
    double r2 = 1.1;
    double c6 = 1.2;
    double c12 = 2.3;
    pele::lj_interaction ljint(c6, c12);
    double h = 0, g = 0;
    double e = ljint.energy_gradient_hessian(r2, &g, &h, 1, 2);
    ASSERT_NEAR(e, 0.39671227804179443, 1e-10);
    ASSERT_NEAR(g, 9.2454671845389917, 1e-10);
    ASSERT_NEAR(h, 149.6972546707778, 1e-10);
}


class LJTest :  public PotentialTest
{
public:

    double c6, c12;
    size_t natoms;

    virtual void setup_potential(){
        pot = std::shared_ptr<pele::BasePotential> (new pele::LJ(c6, c12));
    }

    virtual void SetUp(){
        natoms = 3;
        c6 = 1.2;
        c12 = 2.3;
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
        etrue = -0.10512486948441067;

        setup_potential();
    }


};

TEST_F(LJTest, Energy_Works){
    test_energy();
}

TEST_F(LJTest, EnergyGradient_AgreesWithNumerical){
    test_energy_gradient();
}

TEST_F(LJTest, EnergyGradientHessian_AgreesWithNumerical){
    test_energy_gradient_hessian();
}


/*
 * LJCut
 */


/*
 * LJNeighborList tests
 */

class LJNeighborListTest :  public LJTest
{
    virtual void setup_potential(){
        pele::Array<size_t> ilist(size_t(natoms*(natoms - 1)));
        size_t k = 0;
        for (size_t i=0; i < natoms; ++i){
            for (size_t j=0; j<i; ++j){
                ilist[k] = i;
                ilist[k+1] = j;
                k += 2;
            }
        }
        pot = std::shared_ptr<pele::BasePotential> (new pele::LJNeighborList(ilist, c6, c12));
    }
};

TEST_F(LJNeighborListTest, Energy_Works){
    test_energy();
}

TEST_F(LJNeighborListTest, EnergyGradient_AgreesWithNumerical){
    test_energy_gradient();
}

TEST_F(LJNeighborListTest, EnergyGradientHessian_AgreesWithNumerical){
    test_energy_gradient_hessian();
}


class LJCutTest :  public PotentialTest
{
public:
    double c6, c12, rcut;
    size_t natoms;

    virtual void setup_potential(){
        pot = std::shared_ptr<pele::BasePotential> (new pele::LJCut(c6, c12, rcut));
    }

    virtual void SetUp(){
        c6 = 1.2;
        c12 = 2.3;
        rcut = 2.5;
        natoms = 3;
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

        etrue = -0.089557709975460198;

        setup_potential();
    }
};

TEST_F(LJCutTest, Energy_Works){
    test_energy();
}

TEST_F(LJCutTest, EnergyGradient_AgreesWithNumerical){
    test_energy_gradient();
}

TEST_F(LJCutTest, EnergyGradientHessian_AgreesWithNumerical){
    test_energy_gradient_hessian();
}

/*
 * LJCutAtomList tests
 */

class LJCutAtomListTest :  public LJCutTest
{
public:

    virtual void setup_potential(){
        pele::Array<size_t> atoms(natoms);
        for (size_t i =0; i<atoms.size(); ++i){
            atoms[i] = i;
        }
        pot = std::shared_ptr<pele::BasePotential> (new pele::LJCutAtomList(
                c6, c12, rcut, atoms
                ));
    }

};


TEST_F(LJCutAtomListTest, Energy_Works){
    test_energy();
}

TEST_F(LJCutAtomListTest, EnergyGradient_AgreesWithNumerical){
    test_energy_gradient();
}

TEST_F(LJCutAtomListTest, EnergyGradientHessian_AgreesWithNumerical){
    test_energy_gradient_hessian();
}

class LJCutPeriodicAtomListTest :  public LJCutTest
{
public:

    virtual void setup_potential(){
        //make the box huge, soe the energy is not affected
        //but displace some of the atoms, to check that the distances are being
        //computed correctly
        pele::Array<double> boxvec(3);
        boxvec[0] = 110.;
        boxvec[1] = 121.;
        boxvec[2] = 132.;
        x[0] += 2 * boxvec[0];
        x[4] += 1 * boxvec[1];
        x[5] += 3 * boxvec[2];
        pele::Array<size_t> atoms(natoms);
        for (size_t i =0; i<atoms.size(); ++i){
            atoms[i] = i;
        }
        pot = std::shared_ptr<pele::BasePotential> (new pele::LJCutPeriodicAtomList(
                c6, c12, rcut, boxvec, atoms
                ));
    }

};


TEST_F(LJCutPeriodicAtomListTest, Energy_Works){
    test_energy();
}

TEST_F(LJCutPeriodicAtomListTest, EnergyGradient_AgreesWithNumerical){
    test_energy_gradient();
}

TEST_F(LJCutPeriodicAtomListTest, EnergyGradientHessian_AgreesWithNumerical){
    test_energy_gradient_hessian();
}
