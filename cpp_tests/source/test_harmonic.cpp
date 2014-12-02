#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "pele/harmonic.h"

#include "test_utils.hpp"

using pele::Array;

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
