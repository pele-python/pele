#include <memory>
#include "test_utils.hpp"
#include "pele/lj.h"
#include "pele/lj_cut.h"
#include "pele/atomlist_potential.h"
#include "pele/distance.h"
#include "pele/combine_potentials.h"

class BLJCutTest :  public PotentialTest
{
public:

    double c6, c12, rcut;
    size_t natoms;
    size_t ntypeA;
    pele::Array<size_t> atomsA, atomsB;

    virtual void setup_potential(){
        auto compot = new pele::CombinedPotential();
        compot->add_potential(std::make_shared<pele::LJCutAtomList>(c6, c12, rcut, atomsA));
        compot->add_potential(std::make_shared<pele::LJCutAtomList>(c6, c12, rcut, atomsA, atomsB));
        compot->add_potential(std::make_shared<pele::LJCutAtomList>(c6, c12, rcut, atomsB));

        pot = std::shared_ptr<pele::BasePotential> (compot);

    }

    virtual void SetUp(){
        c6 = 1.2;
        c12 = 2.3;
        rcut = 2.5;
        natoms = 4;
        ntypeA = 2;
        atomsA = pele::Array<size_t>(ntypeA);
        for (size_t i =0; i<atomsA.size(); ++i){
            atomsA[i] = i;
        }
        size_t ntypeB = natoms - ntypeA;
        atomsB = pele::Array<size_t>(ntypeB);
        for (size_t i = 0; i<ntypeB; ++i){
            atomsB[i] = i + ntypeA;
        }

        x = pele::Array<double>(3*natoms);
        x[0]  = 0.1;
        x[1]  = 0.2;
        x[2]  = 0.3;
        x[3]  = 0.44;
        x[4]  = 0.55;
        x[5]  = 1.66;

        x[6] = 0;
        x[7] = 0;
        x[8] = -3.;

        x[9]  = 1.38;
        x[10]  = 0.55;
        x[11]  = 1.66;
        x[9]  = .4;
        x[10]  = 1.57;
        x[11]  = 1.66;

        etrue = 0.65332411282951708;

        setup_potential();
    }



};


TEST_F(BLJCutTest, Energy_Works){
    test_energy();
}

TEST_F(BLJCutTest, EnergyGradient_AgreesWithNumerical){
    test_energy_gradient();
}

TEST_F(BLJCutTest, EnergyGradientHessian_AgreesWithNumerical){
    test_energy_gradient_hessian();
}

