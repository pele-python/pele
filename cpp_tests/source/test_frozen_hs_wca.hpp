#include "pele/array.h"
#include "pele/frozen_atoms.h"
#include "pele/hs_wca.h"

#include <iostream>
#include <stdexcept>
#include <gtest/gtest.h>

using pele::Array;
using pele::HS_WCAFrozen;
using pele::HS_WCAPeriodicFrozen;
using pele::HS_WCA;
using pele::HS_WCAPeriodic;

class FrozenHS_WCATest: public ::testing::Test{
public:
	double eps, sca, etrue;
	Array<double> radii, radii_small, radii_large;
	Array<double> x, y;
	Array<long> frozen_dof;
	double* boxvec;
	virtual void SetUp(){
		eps=1.0;
		sca=1.2;
		radii.resize(3);
		radii[0] = 1.0;
		radii[1] = 1.1;
		radii[2] = 0.9;
		radii_small.resize(3);
		radii_large.resize(3);
		for (unsigned int i(0); i < radii.size(); ++i){
			radii_small[i] = radii[i]/3.;
			radii_large[i] = radii[i]*1.3;
		}
		boxvec = new double[3];
		boxvec[0] = 5;
		boxvec[1] = 6;
		boxvec[2] = 7;
		x.resize(9);
		x[0] = 0.1;
		x[1] = 0.2;
		x[2] = 0.3;
		x[3] = 0.44;
		x[4] = 0.55;
		x[5] = 1.66;
		x[6] = 0.88;
		x[7] = 1.1;
		x[8] = 3.32;
		frozen_dof.resize(3);
		frozen_dof[0] = 0;
		frozen_dof[1] = 3;
		frozen_dof[2] = 4;
	}
	virtual void TearDown() {
		delete[] boxvec;
	}
};

TEST_F(FrozenHS_WCATest, TestEnergy_Correct){
	HS_WCA pot_nofreeze(eps,sca,radii);
	HS_WCAFrozen pot(eps, sca, radii, x, frozen_dof);
	auto xred = pot._coords_converter.get_reduced_coords(x);
	double e = pot.get_energy(xred);
	double etrue =  pot_nofreeze.get_energy(x);
	//std::cout << "etrue: " << etrue << std::endl;
	EXPECT_NEAR(e, etrue, 1e-10);
}

TEST_F(FrozenHS_WCATest, TestEnergy_Correct_Periodic){
	HS_WCAPeriodic pot_nofreeze(eps, sca, radii, boxvec);
	HS_WCAPeriodicFrozen pot(eps, sca, radii, boxvec, x, frozen_dof);
	auto xred = pot._coords_converter.get_reduced_coords(x);
	double e = pot.get_energy(xred);
	double etrue =  pot_nofreeze.get_energy(x);
	EXPECT_NEAR(e, etrue, 1e-10);
}

TEST_F(FrozenHS_WCATest, TestRepulsive_Correct){
	HS_WCA pot_nofreeze(eps,sca,radii);
	HS_WCAFrozen pot(eps, sca, radii, x, frozen_dof);
	HS_WCAFrozen pot_small(eps, sca, radii_small, x, frozen_dof);
	HS_WCAFrozen pot_large(eps, sca, radii_large, x, frozen_dof);
	auto xred = pot._coords_converter.get_reduced_coords(x);
	double e = pot.get_energy(xred);
	double e_small = pot_small.get_energy(xred);
	double e_large = pot_large.get_energy(xred);
	double etrue =  pot_nofreeze.get_energy(x);
	EXPECT_NEAR(e, etrue, 1e-10);
	EXPECT_TRUE(e_small<=e);
	EXPECT_TRUE(e<=e_large);
}

TEST_F(FrozenHS_WCATest, TestRepulsive_Correct_Periodic){
	HS_WCAPeriodic pot_nofreeze(eps, sca, radii, boxvec);
	HS_WCAPeriodicFrozen pot(eps, sca, radii, boxvec, x, frozen_dof);
	HS_WCAPeriodicFrozen pot_small(eps, sca, radii_small, boxvec, x, frozen_dof);
	HS_WCAPeriodicFrozen pot_large(eps, sca, radii_large, boxvec, x, frozen_dof);
	auto xred = pot._coords_converter.get_reduced_coords(x);
	double e = pot.get_energy(xred);
	double e_small = pot_small.get_energy(xred);
	double e_large = pot_large.get_energy(xred);
	double etrue =  pot_nofreeze.get_energy(x);
	EXPECT_NEAR(e, etrue, 1e-10);
	EXPECT_TRUE(e_small<=e);
	EXPECT_TRUE(e<=e_large);
}

TEST_F(FrozenHS_WCATest, TestEnergyGradient_Correct){
	HS_WCA pot_nofreeze(eps,sca,radii);
	HS_WCAFrozen pot(eps, sca, radii, x, frozen_dof);
	auto xred = pot._coords_converter.get_reduced_coords(x);
	Array<double> gred(xred.size()), gtrue(x.size());
	double e = pot.get_energy_gradient(xred,gred);
	double etrue = pot_nofreeze.get_energy_gradient(x,gtrue);
	EXPECT_NEAR(e, etrue, 1e-10);
	auto gtrue_red = pot._coords_converter.get_reduced_coords(gtrue);
	for (int i=0; i<gred.size(); ++i){
		EXPECT_NEAR(gtrue_red[i], gred[i], 1e-10);
	}
}

TEST_F(FrozenHS_WCATest, TestEnergyGradient_Correct_Periodic){
	HS_WCAPeriodic pot_nofreeze(eps, sca, radii, boxvec);
	HS_WCAPeriodicFrozen pot(eps, sca, radii, boxvec, x, frozen_dof);
	auto xred = pot._coords_converter.get_reduced_coords(x);
	Array<double> gred(xred.size()), gtrue(x.size());
	double e = pot.get_energy_gradient(xred,gred);
	double etrue = pot_nofreeze.get_energy_gradient(x,gtrue);
	EXPECT_NEAR(e, etrue, 1e-10);
	auto gtrue_red = pot._coords_converter.get_reduced_coords(gtrue);
	for (int i=0; i<gred.size(); ++i){
		EXPECT_NEAR(gtrue_red[i], gred[i], 1e-10);
	}
}

TEST_F(FrozenHS_WCATest, TestNumericalGradient_Correct){
	HS_WCA pot_nofreeze(eps,sca,radii);
	HS_WCAFrozen pot(eps, sca, radii, x, frozen_dof);
	auto xred = pot._coords_converter.get_reduced_coords(x);
	Array<double> gred(xred.size()), gtrue(x.size());
	pot.numerical_gradient(xred,gred);
	pot_nofreeze.numerical_gradient(x,gtrue);
	auto gtrue_red = pot._coords_converter.get_reduced_coords(gtrue);
	for (int i=0; i<gred.size(); ++i){
		EXPECT_NEAR(gtrue_red[i], gred[i], 1e-10);
	}
}

TEST_F(FrozenHS_WCATest, TestNumericalGradient_Correct_Periodic){
	HS_WCAPeriodic pot_nofreeze(eps, sca, radii, boxvec);
	HS_WCAPeriodicFrozen pot(eps, sca, radii, boxvec, x, frozen_dof);
	auto xred = pot._coords_converter.get_reduced_coords(x);
	Array<double> gred(xred.size()), gtrue(x.size());
	pot.numerical_gradient(xred,gred);
	pot_nofreeze.numerical_gradient(x,gtrue);
	auto gtrue_red = pot._coords_converter.get_reduced_coords(gtrue);
	for (int i=0; i<gred.size(); ++i){
		EXPECT_NEAR(gtrue_red[i], gred[i], 1e-10);
	}
}
