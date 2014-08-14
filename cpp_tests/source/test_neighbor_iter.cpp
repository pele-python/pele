#include "pele/array.h"
#include "pele/inversepower.h"
#include "pele/neighbor_iterator.h"

#include <iostream>
#include <stdexcept>
#include <gtest/gtest.h>

using pele::Array;
using pele::InversePowerPeriodic;
using pele::InversePower_interaction;

/*
 * HS_WCA tests
 */

class CellIterTest :  public ::testing::Test
{
public:
    double pow, eps, etrue, rcut;
    Array<double> x, g, gnum, radii, boxvec;
    virtual void SetUp(){
    	pow = 2.5;
    	eps = 1;
    	rcut = 4;
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
        boxvec = Array<double>(3,rcut);
        double f = 1.;
        radii[0] = .91 * f;
        radii[1] = 1.1 * f;
        radii[2] = 1.13 * f;
        etrue = 0.03493116137645523;
        g = Array<double>(x.size());
        gnum = Array<double>(x.size());
    }
};

//test number of distinguishable pairs
TEST_F(CellIterTest, Number_of_neighbors){
    pele::CellIter<> cell = pele::CellIter<>(x, boxvec, rcut);
    pele::CellIter<>::const_iterator it;
    size_t count = 0;
    for(it = cell.begin(); it != cell.end(); ++it){
        ++count;
    }
    ASSERT_EQ(3u, count);
    ASSERT_EQ(3u, static_cast<unsigned int>(cell.end() - cell.begin()));
    ASSERT_EQ(3u, cell.get_nr_unique_pairs());
}

/*TEST_F(CellIterTest, Energy_Works){
    InversePowerPeriodic<3> pot(pow, eps, radii, boxvec);
    double e = pot.get_energy(x);

    double ecell =0;
    for(int i=0;i<3;++i){
        ecell +=
    }
    ASSERT_NEAR(e, etrue, 1e-10);
}

TEST_F(CellIterTest, EnergyGradient_AgreesWithNumerical){
	InversePowerPeriodic<3> pot(pow, eps, radii, boxvec);
	double e = pot.get_energy_gradient(x, g);
	std::cout<<"energy"<<e<<std::endl;
    double ecomp = pot.get_energy(x);
    ASSERT_NEAR(e, ecomp, 1e-10);
    pot.numerical_gradient(x, gnum, 1e-6);
    for (size_t k=0; k<6; ++k){
        ASSERT_NEAR(g[k], gnum[k], 1e-6);
    }
}

TEST_F(CellIterTest, EnergyGradientHessian_AgreesWithNumerical){
	InversePowerPeriodic<3> pot(pow, eps, radii, boxvec);
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
}*/
