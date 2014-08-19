#include "pele/array.h"
#include "pele/inversepower.h"
#include "pele/hs_wca.h"
#include "pele/neighbor_iterator.h"
#include "pele/modified_fire.h"
#include "pele/distance.h"

#include <iostream>
#include <stdexcept>
#include <gtest/gtest.h>
#include <random>
#include <ctime>

using pele::Array;
using pele::InversePowerPeriodic;
using pele::InversePower_interaction;

class CellIterTest : public ::testing::Test {
public:
    double pow, eps, etrue, rcut, sca;
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
        boxvec = Array<double>(3, rcut);
        double f = 1.;
        radii[0] = .91 * f;
        radii[1] = 1.1 * f;
        radii[2] = 1.13 * f;
        etrue = 0.03493116137645523;
        g = Array<double>(x.size());
        gnum = Array<double>(x.size());
        sca = 1.2;
    }
};

//test number of distinguishable pairs
TEST_F(CellIterTest, Number_of_neighbors){
    pele::CellIter<> cell = pele::CellIter<>(x, boxvec, rcut);
    pele::CellIter<> cell2 = pele::CellIter<>(x, boxvec, rcut, 1);
    pele::CellIter<> cell3 = pele::CellIter<>(x, boxvec, rcut, 4.2);
    pele::CellIter<> cell4 = pele::CellIter<>(x, boxvec, rcut, 5);
    size_t count = 0;
    size_t count2 = 0;
    size_t count3 = 0;
    size_t count4 = 0;
    pele::CellIter<>::const_iterator it;
    for (it = cell.begin(); it != cell.end(); ++it, ++count);
    for (it = cell2.begin(); it != cell2.end(); ++it, ++count2);
    for (it = cell3.begin(); it != cell3.end(); ++it, ++count3);
    for (it = cell4.begin(); it != cell4.end(); ++it, ++count4);
    ASSERT_EQ(3u, count);
    ASSERT_EQ(count, static_cast<unsigned int>(cell.end() - cell.begin()));
    ASSERT_EQ(count, cell.get_nr_unique_pairs());
    ASSERT_EQ(count, count2);
    ASSERT_EQ(count, static_cast<unsigned int>(cell2.end() - cell2.begin()));
    ASSERT_EQ(count, cell2.get_nr_unique_pairs());
    ASSERT_EQ(count, count3);
    ASSERT_EQ(count, static_cast<unsigned int>(cell3.end() - cell3.begin()));
    ASSERT_EQ(count, cell3.get_nr_unique_pairs());
    ASSERT_EQ(count, count4);
    ASSERT_EQ(count, static_cast<unsigned int>(cell4.end() - cell4.begin()));
    ASSERT_EQ(count, cell4.get_nr_unique_pairs());
}

TEST_F(CellIterTest, Energy_Works){
    pele::InversePowerPeriodicCellLists<3> pot_cell(pow, eps, radii, boxvec, x, rcut, 1.0);
    pele::InversePowerPeriodicCellLists<3> pot_cell2(pow, eps, radii, boxvec, x, rcut, 2.0);
    pele::InversePowerPeriodicCellLists<3> pot_cell3(pow, eps, radii, boxvec, x, rcut, 3.0);
    pele::InversePowerPeriodicCellLists<3> pot_cell4(pow, eps, radii, boxvec, x, rcut, 4.0);
    pele::InversePowerPeriodic<3> pot(pow, eps, radii, boxvec);
    const double ecell = pot_cell.get_energy(x);
    const double ecell2 = pot_cell2.get_energy(x);
    const double ecell3 = pot_cell3.get_energy(x);
    const double ecell4 = pot_cell4.get_energy(x);
    const double etrue = pot.get_energy(x);
    ASSERT_NEAR(ecell, etrue, 1e-10);
    ASSERT_NEAR(ecell2, etrue, 1e-10);
    ASSERT_NEAR(ecell3, etrue, 1e-10);
    ASSERT_NEAR(ecell4, etrue, 1e-10);
}

TEST_F(CellIterTest, EnergyGradient_AgreesWithNumerical){
    pele::InversePowerPeriodic<3> pot_no_cells(pow, eps, radii, boxvec);
    const double etrue = pot_no_cells.get_energy(x);
    const size_t N = 4;
    std::vector<std::shared_ptr<pele::InversePowerPeriodicCellLists<3> > > pot;
    for (size_t i = 0; i < N; ++i) {
        pot.push_back(std::make_shared<pele::InversePowerPeriodicCellLists<3> >(
                pow, eps, radii, boxvec, x, rcut, 1 + i));
    }
    pot.swap(pot);
    std::vector<double> e(N, 0);
    std::vector<double> ecomp(N, 0);
    for (size_t i = 0; i < N; ++i) {
        e.at(i) = pot.at(i)->get_energy_gradient(x, g);
        ecomp.at(i) = pot.at(i)->get_energy(x);
        pot.at(i)->numerical_gradient(x, gnum, 1e-6);
        for (size_t k = 0; k < 6; ++k) {
            ASSERT_NEAR(g[k], gnum[k], 1e-6);
        }
    }
    for (size_t i = 0; i < N; ++i) {
        ASSERT_NEAR(e.at(i), ecomp.at(i), 1e-10);
        ASSERT_NEAR(e.at(i), etrue, 1e-10);
    }
}

TEST_F(CellIterTest, EnergyGradientHessian_AgreesWithNumerical){
    pele::InversePowerPeriodic<3> pot_no_cells(pow, eps, radii, boxvec);
    const double etrue = pot_no_cells.get_energy(x);
    Array<double> g_no_cells(x.size()) ;
    Array<double> h_no_cells(x.size() * x.size());
    pot_no_cells.get_energy_gradient_hessian(x, g_no_cells, h_no_cells);
    for (size_t i = 0; i < 4; ++i) {
        pele::InversePowerPeriodicCellLists<3> pot(pow, eps, radii, boxvec, x, rcut, 1.0 + i);
        Array<double> h(x.size() * x.size());
        Array<double> hnum(h.size());
        const double e = pot.get_energy_gradient_hessian(x, g, h);
        const double ecomp = pot.get_energy(x);
        pot.numerical_gradient(x, gnum);
        pot.numerical_hessian(x, hnum);
        EXPECT_NEAR(e, ecomp, 1e-10);
        EXPECT_NEAR(etrue, ecomp, 1e-10);
        for (size_t i = 0; i < g.size(); ++i) {
            ASSERT_NEAR(g[i], gnum[i], 1e-10);
            ASSERT_NEAR(g[i], g_no_cells[i], 1e-10);
        }
        for (size_t i = 0; i < h.size(); ++i) {
            ASSERT_NEAR(h[i], hnum[i], 1e-10);
            ASSERT_NEAR(h[i], h_no_cells[i], 1e-10);
        }
    }
}

TEST_F(CellIterTest, HS_WCAEnergy_Works){
    pele::HS_WCAPeriodicCellLists<3> pot_cell(eps, sca, radii, boxvec, x, rcut, 1);
    pele::HS_WCAPeriodicCellLists<3> pot_cell2(eps, sca, radii, boxvec, x, rcut, 2);
    pele::HS_WCAPeriodicCellLists<3> pot_cell3(eps, sca, radii, boxvec, x, rcut, 3);
    pele::HS_WCAPeriodicCellLists<3> pot_cell4(eps, sca, radii, boxvec, x, rcut, 4);
    pele::HS_WCAPeriodic<3> pot_no_cells(eps, sca, radii, boxvec);
    const double ecell = pot_cell.get_energy(x);
    const double ecell2 = pot_cell2.get_energy(x);
    const double ecell3 = pot_cell3.get_energy(x);
    const double ecell4 = pot_cell4.get_energy(x);
    const double etrue = pot_no_cells.get_energy(x);
    ASSERT_NEAR(ecell, etrue, 1e-10);
    ASSERT_NEAR(ecell2, etrue, 1e-10);
    ASSERT_NEAR(ecell3, etrue, 1e-10);
    ASSERT_NEAR(ecell4, etrue, 1e-10);
}

class CellIterTestMoreHS_WCA : public ::testing::Test {
public:
    double pow;
    size_t seed;
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> distribution;
    size_t nparticles;
    size_t ndim;
    size_t ndof;
    double eps;
    double rcut;
    double sca;
    Array<double> x;
    Array<double> g;
    Array<double> gnum;
    Array<double> radii;
    Array<double> boxvec;
    virtual void SetUp(){
        pow = 2.5;
        seed = 42;
        generator = std::mt19937_64(seed);
        distribution = std::uniform_real_distribution<double>(0, 0.05);
        nparticles = 20;
        ndim = 3;
        ndof = nparticles * ndim;
        eps = 1;
        rcut = 4;
        x = Array<double>(ndof);
        for (size_t k = 0; k < ndof; ++k) {
            x[k] = double(k + 1) / double(10) + (k / nparticles) * (1 + distribution(generator));
        }
        radii = Array<double>(nparticles);
        boxvec = Array<double>(ndim, 9);
        for (size_t i = 0; i < nparticles; ++i) {
            radii[i] = (0.1 + distribution(generator));
        }
        g = Array<double>(x.size());
        gnum = Array<double>(x.size());
        sca = 1.2;
    }
};

TEST_F(CellIterTestMoreHS_WCA, Number_of_neighbors){
    pele::CellIter<> cell = pele::CellIter<>(x, boxvec, boxvec[0]);
    pele::CellIter<> cell2 = pele::CellIter<>(x, boxvec, boxvec[0], 1);
    pele::CellIter<> cell3 = pele::CellIter<>(x, boxvec, boxvec[0], 3);
    pele::CellIter<> cell4 = pele::CellIter<>(x, boxvec, boxvec[0], 5);
    size_t count = 0;
    size_t count2 = 0;
    size_t count3 = 0;
    size_t count4 = 0;
    pele::CellIter<>::const_iterator it;
    for (it = cell.begin(); it != cell.end(); ++it, ++count);
    for (it = cell2.begin(); it != cell2.end(); ++it, ++count2);
    for (it = cell3.begin(); it != cell3.end(); ++it, ++count3);
    for (it = cell4.begin(); it != cell4.end(); ++it, ++count4);
    ASSERT_EQ(nparticles * (nparticles - 1) / 2, count);
    ASSERT_EQ(count, static_cast<unsigned int>(cell.end() - cell.begin()));
    ASSERT_EQ(count, cell.get_nr_unique_pairs());
    ASSERT_EQ(count, count2);
    ASSERT_EQ(count, static_cast<unsigned int>(cell2.end() - cell2.begin()));
    ASSERT_EQ(count, cell2.get_nr_unique_pairs());
    ASSERT_EQ(count, count3);
    ASSERT_EQ(count, static_cast<unsigned int>(cell3.end() - cell3.begin()));
    ASSERT_EQ(count, cell3.get_nr_unique_pairs());
    ASSERT_EQ(count, count4);
    ASSERT_EQ(count, static_cast<unsigned int>(cell4.end() - cell4.begin()));
    ASSERT_EQ(count, cell4.get_nr_unique_pairs());
}

TEST_F(CellIterTestMoreHS_WCA, EnergyMoreParticles_Works){
    pele::InversePowerPeriodicCellLists<3> pot_cell(pow, eps, radii, boxvec, x, rcut, 1);
    pele::InversePowerPeriodicCellLists<3> pot_cell2(pow, eps, radii, boxvec, x, rcut, 2);
    pele::InversePowerPeriodicCellLists<3> pot_cell3(pow, eps, radii, boxvec, x, rcut, 3);
    pele::InversePowerPeriodic<3> pot(pow, eps, radii, boxvec);
    const double ecell = pot_cell.get_energy(x);
    const double ecell2 = pot_cell2.get_energy(x);
    const double ecell3 = pot_cell3.get_energy(x);
    const double etrue = pot.get_energy(x);
    EXPECT_DOUBLE_EQ(ecell, etrue);
    EXPECT_DOUBLE_EQ(ecell2, etrue);
    EXPECT_DOUBLE_EQ(ecell3, etrue);
}

TEST_F(CellIterTestMoreHS_WCA, HSWCAEnergy_Works) {
    pele::HS_WCAPeriodic<3> pot_no_cells(eps, sca, radii, boxvec);
    const double e_no_cells = pot_no_cells.get_energy(x);
    for (size_t factor = 1; factor < 4; ++factor) {
        pele::HS_WCAPeriodicCellLists<3> pot_cellA(eps, sca, radii, boxvec, x, rcut, factor);
        pele::HS_WCAPeriodicCellLists<3> pot_cellB(eps, sca, radii, boxvec, x, rcut, factor + 0.2);
        const double e_cellA = pot_cellA.get_energy(x);
        const double e_cellB = pot_cellB.get_energy(x);
        EXPECT_DOUBLE_EQ(e_no_cells, e_cellA);
        EXPECT_DOUBLE_EQ(e_no_cells, e_cellB);
    }
}

TEST_F(CellIterTestMoreHS_WCA, HSWCAEnergyGradient_Works) {
    pele::HS_WCAPeriodic<3> pot_no_cells(eps, sca, radii, boxvec);
    const double e_no_cells = pot_no_cells.get_energy(x);
    for (size_t factor = 1; factor < 4; ++factor) {
        pele::HS_WCAPeriodicCellLists<3> pot_cellA(eps, sca, radii, boxvec, x, rcut, factor);
        pele::HS_WCAPeriodicCellLists<3> pot_cellB(eps, sca, radii, boxvec, x, rcut, factor + 0.2);
        pele::Array<double> g_no_cells(x.size());
        pele::Array<double> g_cellA(x.size());
        pele::Array<double> g_cellB(x.size());
        const double eg_no_cells = pot_no_cells.get_energy_gradient(x, g_no_cells);
        const double eg_cellA = pot_cellA.get_energy_gradient(x, g_cellA);
        const double eg_cellB = pot_cellB.get_energy_gradient(x, g_cellB);
        EXPECT_DOUBLE_EQ(e_no_cells, eg_no_cells);
        EXPECT_DOUBLE_EQ(e_no_cells, eg_cellA);
        EXPECT_DOUBLE_EQ(e_no_cells, eg_cellB);
        for (size_t i = 0; i < g_no_cells.size(); ++i) {
            EXPECT_DOUBLE_EQ(g_no_cells[i], g_cellA[i]);
            EXPECT_DOUBLE_EQ(g_no_cells[i], g_cellB[i]);
        }
    }
}

TEST_F(CellIterTestMoreHS_WCA, HSWCAEnergyGradientHessian_Works) {
    pele::HS_WCAPeriodic<3> pot_no_cells(eps, sca, radii, boxvec);
    const double e_no_cells = pot_no_cells.get_energy(x);
    for (size_t factor = 1; factor < 4; ++factor) {
        pele::HS_WCAPeriodicCellLists<3> pot_cellA(eps, sca, radii, boxvec, x, rcut, factor);
        pele::HS_WCAPeriodicCellLists<3> pot_cellB(eps, sca, radii, boxvec, x, rcut, factor + 0.2);
        pele::Array<double> g_no_cells(x.size());
        pele::Array<double> g_cellA(x.size());
        pele::Array<double> g_cellB(x.size());
        pele::Array<double> h_no_cells(x.size() * x.size());
        pele::Array<double> h_cellA(h_no_cells.size());
        pele::Array<double> h_cellB(h_no_cells.size());
        const double egh_no_cells = pot_no_cells.get_energy_gradient_hessian(x, g_no_cells, h_no_cells);
        const double egh_cellA = pot_cellA.get_energy_gradient_hessian(x, g_cellA, h_cellA);
        const double egh_cellB = pot_cellB.get_energy_gradient_hessian(x, g_cellB, h_cellB);
        EXPECT_DOUBLE_EQ(e_no_cells, egh_no_cells);
        EXPECT_DOUBLE_EQ(e_no_cells, egh_cellA);
        EXPECT_DOUBLE_EQ(e_no_cells, egh_cellB);
        for (size_t i = 0; i < g_no_cells.size(); ++i) {
            EXPECT_DOUBLE_EQ(g_no_cells[i], g_cellA[i]);
            EXPECT_DOUBLE_EQ(g_no_cells[i], g_cellB[i]);
        }
        for (size_t i = 0; i < h_no_cells.size(); ++i) {
            EXPECT_DOUBLE_EQ(h_no_cells[i], h_cellA[i]);
            EXPECT_DOUBLE_EQ(h_no_cells[i], h_cellB[i]);
        }
    }
}

TEST_F(CellIterTestMoreHS_WCA, HSWCAMinimzation_Works) {
    auto pot_no_cells = std::make_shared<pele::HS_WCAPeriodic<3> >(eps, sca, radii, boxvec);
    auto pot_cells = std::make_shared<pele::HS_WCAPeriodicCellLists<3> >(eps, sca, radii, boxvec, x, rcut, 1);
    pele::MODIFIED_FIRE opt_no_cells(pot_no_cells, x, .1, 1, 1);
    pele::MODIFIED_FIRE opt_cells(pot_cells, x, .1, 1, 1);
    opt_no_cells.run();
    opt_cells.run();
    auto x_opt_no_cells = opt_no_cells.get_x();
    auto x_opt_cells = opt_no_cells.get_x();
    const auto e_opt_no_cells = pot_no_cells->get_energy(x_opt_no_cells);
    const auto e_opt_cells = pot_cells->get_energy(x_opt_cells);
    EXPECT_DOUBLE_EQ(e_opt_no_cells, e_opt_cells);
    for (size_t i = 0; i < x_opt_no_cells.size(); ++i) {
        EXPECT_DOUBLE_EQ(x_opt_no_cells[i], x_opt_cells[i]);
    }
}

TEST_F(CellIterTestMoreHS_WCA, HSWCAMinimzationPerformance_Changes) {
    auto pot_no_cells = std::make_shared<pele::HS_WCAPeriodic<3> >(eps, sca, radii, boxvec);
    auto pot_cells = std::make_shared<pele::HS_WCAPeriodicCellLists<3> >(eps, sca, radii, boxvec, x, 0.4, 0.3);
    pele::MODIFIED_FIRE opt_no_cells(pot_no_cells, x, .1, 1, 1);
    pele::MODIFIED_FIRE opt_cells(pot_cells, x, .1, 1, 1);
    const size_t repetitions = 50;
    const auto start_time = clock();
    for (size_t i = 0; i < repetitions; ++i) {
        opt_no_cells.reset(x);
        opt_no_cells.run();
    }
    const auto intermediate_1 = clock();
    for (size_t i = 0; i < repetitions; ++i) {
        opt_cells.reset(x);
        opt_cells.run();
    }
    const auto intermediate_2 = clock();
    auto x_opt_no_cells = opt_no_cells.get_x();
    auto x_opt_cells = opt_no_cells.get_x();
    const auto e_opt_no_cells = pot_no_cells->get_energy(x_opt_no_cells);
    const auto e_opt_cells = pot_cells->get_energy(x_opt_cells);
    EXPECT_DOUBLE_EQ(e_opt_no_cells, e_opt_cells);
    for (size_t i = 0; i < x_opt_no_cells.size(); ++i) {
        EXPECT_DOUBLE_EQ(x_opt_no_cells[i], x_opt_cells[i]);
    }
    const auto time_no_cells = intermediate_1 - start_time;
    const auto time_cells = intermediate_2 - intermediate_1;
    EXPECT_LE(time_cells, time_no_cells);
}

class CellIterTestMoreHS_WCA2D : public ::testing::Test {
public:
    size_t seed;
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> distribution;
    size_t nparticles;
    size_t ndim;
    size_t ndof;
    double eps;
    double rcut;
    double sca;
    Array<double> x;
    Array<double> g;
    Array<double> gnum;
    Array<double> radii;
    Array<double> boxvec;
    virtual void SetUp(){
        seed = 42;
        generator = std::mt19937_64(seed);
        distribution = std::uniform_real_distribution<double>(0, 0.05);
        nparticles = 20;
        ndim = 2;
        ndof = nparticles * ndim;
        eps = 1;
        rcut = 4;
        x = Array<double>(ndof);
        for (size_t k = 0; k < ndof; ++k) {
            x[k] = double(k + 1) / double(10) + (k / nparticles) * (1 + distribution(generator));
            x[k] -= distribution(generator) * 0.5;
        }
        radii = Array<double>(nparticles);
        boxvec = Array<double>(ndim, 9);
        for (size_t i = 0; i < nparticles; ++i) {
            radii[i] = (0.1 + distribution(generator));
        }
        g = Array<double>(x.size());
        gnum = Array<double>(x.size());
        sca = 1.2;
    }
};

TEST_F(CellIterTestMoreHS_WCA2D, Number_of_neighbors){
    auto cell = pele::CellIter<pele::periodic_distance<2> >(x, boxvec, boxvec[0]);
    auto cell2 = pele::CellIter<pele::periodic_distance<2> >(x, boxvec, boxvec[0], 1);
    auto cell3 = pele::CellIter<pele::periodic_distance<2> >(x, boxvec, boxvec[0], 3);
    auto cell4 = pele::CellIter<pele::periodic_distance<2> >(x, boxvec, boxvec[0], 5);
    size_t count = 0;
    size_t count2 = 0;
    size_t count3 = 0;
    size_t count4 = 0;
    pele::CellIter<pele::periodic_distance<2> >::const_iterator it;
    for (it = cell.begin(); it != cell.end(); ++it, ++count);
    for (it = cell2.begin(); it != cell2.end(); ++it, ++count2);
    for (it = cell3.begin(); it != cell3.end(); ++it, ++count3);
    for (it = cell4.begin(); it != cell4.end(); ++it, ++count4);
    ASSERT_EQ(nparticles * (nparticles - 1) / 2, count);
    ASSERT_EQ(count, static_cast<unsigned int>(cell.end() - cell.begin()));
    ASSERT_EQ(count, cell.get_nr_unique_pairs());
    ASSERT_EQ(count, count2);
    ASSERT_EQ(count, static_cast<unsigned int>(cell2.end() - cell2.begin()));
    ASSERT_EQ(count, cell2.get_nr_unique_pairs());
    ASSERT_EQ(count, count3);
    ASSERT_EQ(count, static_cast<unsigned int>(cell3.end() - cell3.begin()));
    ASSERT_EQ(count, cell3.get_nr_unique_pairs());
    ASSERT_EQ(count, count4);
    ASSERT_EQ(count, static_cast<unsigned int>(cell4.end() - cell4.begin()));
    ASSERT_EQ(count, cell4.get_nr_unique_pairs());
}

TEST_F(CellIterTestMoreHS_WCA2D, HSWCAEnergy_Works) {
    pele::HS_WCAPeriodic<2> pot_no_cells(eps, sca, radii, boxvec);
    const double e_no_cells = pot_no_cells.get_energy(x);
    for (size_t factor = 1; factor < 4; ++factor) {
        pele::HS_WCAPeriodicCellLists<2> pot_cellA(eps, sca, radii, boxvec, x, rcut, factor);
        pele::HS_WCAPeriodicCellLists<2> pot_cellB(eps, sca, radii, boxvec, x, rcut, factor + 0.2);
        const double e_cellA = pot_cellA.get_energy(x);
        const double e_cellB = pot_cellB.get_energy(x);
        EXPECT_DOUBLE_EQ(e_no_cells, e_cellA);
        EXPECT_DOUBLE_EQ(e_no_cells, e_cellB);
    }
}

TEST_F(CellIterTestMoreHS_WCA2D, HSWCAEnergyGradient_Works) {
    pele::HS_WCAPeriodic<2> pot_no_cells(eps, sca, radii, boxvec);
    const double e_no_cells = pot_no_cells.get_energy(x);
    for (size_t factor = 1; factor < 4; ++factor) {
        pele::HS_WCAPeriodicCellLists<2> pot_cellA(eps, sca, radii, boxvec, x, rcut, factor);
        pele::HS_WCAPeriodicCellLists<2> pot_cellB(eps, sca, radii, boxvec, x, rcut, factor + 0.2);
        pele::Array<double> g_no_cells(x.size());
        pele::Array<double> g_cellA(x.size());
        pele::Array<double> g_cellB(x.size());
        const double eg_no_cells = pot_no_cells.get_energy_gradient(x, g_no_cells);
        const double eg_cellA = pot_cellA.get_energy_gradient(x, g_cellA);
        const double eg_cellB = pot_cellB.get_energy_gradient(x, g_cellB);
        EXPECT_DOUBLE_EQ(e_no_cells, eg_no_cells);
        EXPECT_DOUBLE_EQ(e_no_cells, eg_cellA);
        EXPECT_DOUBLE_EQ(e_no_cells, eg_cellB);
        for (size_t i = 0; i < g_no_cells.size(); ++i) {
            EXPECT_DOUBLE_EQ(g_no_cells[i], g_cellA[i]);
            EXPECT_DOUBLE_EQ(g_no_cells[i], g_cellB[i]);
        }
    }
}

TEST_F(CellIterTestMoreHS_WCA2D, HSWCAEnergyGradientHessian_Works) {
    pele::HS_WCAPeriodic<2> pot_no_cells(eps, sca, radii, boxvec);
    const double e_no_cells = pot_no_cells.get_energy(x);
    for (size_t factor = 1; factor < 4; ++factor) {
        pele::HS_WCAPeriodicCellLists<2> pot_cellA(eps, sca, radii, boxvec, x, rcut, factor);
        pele::HS_WCAPeriodicCellLists<2> pot_cellB(eps, sca, radii, boxvec, x, rcut, factor + 0.2);
        pele::Array<double> g_no_cells(x.size());
        pele::Array<double> g_cellA(x.size());
        pele::Array<double> g_cellB(x.size());
        pele::Array<double> h_no_cells(x.size() * x.size());
        pele::Array<double> h_cellA(h_no_cells.size());
        pele::Array<double> h_cellB(h_no_cells.size());
        const double egh_no_cells = pot_no_cells.get_energy_gradient_hessian(x, g_no_cells, h_no_cells);
        const double egh_cellA = pot_cellA.get_energy_gradient_hessian(x, g_cellA, h_cellA);
        const double egh_cellB = pot_cellB.get_energy_gradient_hessian(x, g_cellB, h_cellB);
        EXPECT_DOUBLE_EQ(e_no_cells, egh_no_cells);
        EXPECT_DOUBLE_EQ(e_no_cells, egh_cellA);
        EXPECT_DOUBLE_EQ(e_no_cells, egh_cellB);
        for (size_t i = 0; i < g_no_cells.size(); ++i) {
            EXPECT_DOUBLE_EQ(g_no_cells[i], g_cellA[i]);
            EXPECT_DOUBLE_EQ(g_no_cells[i], g_cellB[i]);
        }
        for (size_t i = 0; i < h_no_cells.size(); ++i) {
            EXPECT_DOUBLE_EQ(h_no_cells[i], h_cellA[i]);
            EXPECT_DOUBLE_EQ(h_no_cells[i], h_cellB[i]);
        }
    }
}

TEST_F(CellIterTestMoreHS_WCA2D, HSWCAMinimzation_Works) {
    auto pot_no_cells = std::make_shared<pele::HS_WCAPeriodic<2> >(eps, sca, radii, boxvec);
    auto pot_cells = std::make_shared<pele::HS_WCAPeriodicCellLists<2> >(eps, sca, radii, boxvec, x, rcut, 1);
    pele::MODIFIED_FIRE opt_no_cells(pot_no_cells, x, .1, 1, 1);
    pele::MODIFIED_FIRE opt_cells(pot_cells, x, .1, 1, 1);
    opt_no_cells.run();
    opt_cells.run();
    auto x_opt_no_cells = opt_no_cells.get_x();
    auto x_opt_cells = opt_no_cells.get_x();
    const auto e_opt_no_cells = pot_no_cells->get_energy(x_opt_no_cells);
    const auto e_opt_cells = pot_cells->get_energy(x_opt_cells);
    EXPECT_DOUBLE_EQ(e_opt_no_cells, e_opt_cells);
    for (size_t i = 0; i < x_opt_no_cells.size(); ++i) {
        EXPECT_DOUBLE_EQ(x_opt_no_cells[i], x_opt_cells[i]);
    }
}

TEST_F(CellIterTestMoreHS_WCA2D, HSWCAMinimzationPerformance_Changes) {
    auto pot_no_cells = std::make_shared<pele::HS_WCAPeriodic<2> >(eps, sca, radii, boxvec);
    auto pot_cells = std::make_shared<pele::HS_WCAPeriodicCellLists<2> >(eps, sca, radii, boxvec, x, 1, 0.3);
    pele::MODIFIED_FIRE opt_no_cells(pot_no_cells, x, .1, 1, 1);
    pele::MODIFIED_FIRE opt_cells(pot_cells, x, .1, 1, 1);
    const size_t repetitions = 50;
    const auto start_time = clock();
    for (size_t i = 0; i < repetitions; ++i) {
        opt_no_cells.reset(x);
        opt_no_cells.run();
    }
    const auto intermediate_1 = clock();
    for (size_t i = 0; i < repetitions; ++i) {
        opt_cells.reset(x);
        opt_cells.run();
    }
    const auto intermediate_2 = clock();
    auto x_opt_no_cells = opt_no_cells.get_x();
    auto x_opt_cells = opt_no_cells.get_x();
    const auto e_opt_no_cells = pot_no_cells->get_energy(x_opt_no_cells);
    const auto e_opt_cells = pot_cells->get_energy(x_opt_cells);
    EXPECT_DOUBLE_EQ(e_opt_no_cells, e_opt_cells);
    for (size_t i = 0; i < x_opt_no_cells.size(); ++i) {
        EXPECT_DOUBLE_EQ(x_opt_no_cells[i], x_opt_cells[i]);
    }
    const auto time_no_cells = intermediate_1 - start_time;
    const auto time_cells = intermediate_2 - intermediate_1;
    EXPECT_LE(time_cells, time_no_cells);
}

class CellIterTestMoreHS_WCA2DFrozen : public ::testing::Test {
public:
    size_t seed;
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> distribution;
    size_t nparticles;
    size_t ndim;
    size_t ndof;
    size_t n_frozen_dof;
    std::uniform_int_distribution<size_t> int_distr;
    double eps;
    double rcut;
    double sca;
    Array<double> x;
    Array<double> g;
    Array<double> gnum;
    Array<double> radii;
    Array<double> boxvec;
    Array<size_t> frozen_dof;
    virtual void SetUp(){
        seed = 42;
        generator = std::mt19937_64(seed);
        distribution = std::uniform_real_distribution<double>(0, 0.05);
        nparticles = 20;
        ndim = 2;
        ndof = nparticles * ndim;
        n_frozen_dof = ndof * 0.25;
        int_distr = std::uniform_int_distribution<size_t>(0, ndof - 1);
        frozen_dof = Array<size_t>(n_frozen_dof);
        for (size_t i = 0; i < n_frozen_dof; ++i) {
            frozen_dof[i] = int_distr(generator);
        }
        eps = 1;
        rcut = 4;
        x = Array<double>(ndof);
        for (size_t k = 0; k < ndof; ++k) {
            x[k] = double(k + 1) / double(10) + (k / nparticles) * (1 + distribution(generator));
            x[k] -= distribution(generator) * 0.5;
        }
        radii = Array<double>(nparticles);
        boxvec = Array<double>(ndim, 9);
        for (size_t i = 0; i < nparticles; ++i) {
            radii[i] = (0.1 + distribution(generator));
        }
        g = Array<double>(x.size());
        gnum = Array<double>(x.size());
        sca = 1.2;
    }
};

TEST_F(CellIterTestMoreHS_WCA2DFrozen, HSWCAEnergy_Works) {
    pele::HS_WCAPeriodic<2> pot_no_cells_nofr(eps, sca, radii, boxvec);
    pele::HS_WCAPeriodicFrozen<2> pot_no_cells(eps, sca, radii, boxvec, x, frozen_dof);
    auto xred_no_cells = pot_no_cells.coords_converter.get_reduced_coords(x);
    const double e_no_cells = pot_no_cells.get_energy(xred_no_cells);
    const double e_true = pot_no_cells_nofr.get_energy(x);
    EXPECT_DOUBLE_EQ(e_true, e_no_cells);
    for (size_t factor = 1; factor < 4; ++factor) {
        pele::HS_WCAPeriodicCellListsFrozen<2> pot_cellA(eps, sca, radii, boxvec, x, frozen_dof, rcut, factor);
        pele::HS_WCAPeriodicCellListsFrozen<2> pot_cellB(eps, sca, radii, boxvec, x, frozen_dof, rcut, factor + 0.2);
        auto xred_cellA = pot_cellA.coords_converter.get_reduced_coords(x);
        auto xred_cellB = pot_cellB.coords_converter.get_reduced_coords(x);
        const double e_cellA = pot_cellA.get_energy(xred_cellA);
        const double e_cellB = pot_cellB.get_energy(xred_cellB);
        EXPECT_DOUBLE_EQ(e_no_cells, e_cellA);
        EXPECT_DOUBLE_EQ(e_no_cells, e_cellB);
    }
}





















////////////////////////////
////////////////////////////

/*




TEST_F(CellIterTestMoreHS_WCA2D, HSWCAEnergyGradient_Works) {
    pele::HS_WCAPeriodic<2> pot_no_cells(eps, sca, radii, boxvec);
    const double e_no_cells = pot_no_cells.get_energy(x);
    for (size_t factor = 1; factor < 4; ++factor) {
        pele::HS_WCAPeriodicCellLists<2> pot_cellA(eps, sca, radii, boxvec, x, rcut, factor);
        pele::HS_WCAPeriodicCellLists<2> pot_cellB(eps, sca, radii, boxvec, x, rcut, factor + 0.2);
        pele::Array<double> g_no_cells(x.size());
        pele::Array<double> g_cellA(x.size());
        pele::Array<double> g_cellB(x.size());
        const double eg_no_cells = pot_no_cells.get_energy_gradient(x, g_no_cells);
        const double eg_cellA = pot_cellA.get_energy_gradient(x, g_cellA);
        const double eg_cellB = pot_cellB.get_energy_gradient(x, g_cellB);
        EXPECT_DOUBLE_EQ(e_no_cells, eg_no_cells);
        EXPECT_DOUBLE_EQ(e_no_cells, eg_cellA);
        EXPECT_DOUBLE_EQ(e_no_cells, eg_cellB);
        for (size_t i = 0; i < g_no_cells.size(); ++i) {
            EXPECT_DOUBLE_EQ(g_no_cells[i], g_cellA[i]);
            EXPECT_DOUBLE_EQ(g_no_cells[i], g_cellB[i]);
        }
    }
}

TEST_F(CellIterTestMoreHS_WCA2D, HSWCAEnergyGradientHessian_Works) {
    pele::HS_WCAPeriodic<2> pot_no_cells(eps, sca, radii, boxvec);
    const double e_no_cells = pot_no_cells.get_energy(x);
    for (size_t factor = 1; factor < 4; ++factor) {
        pele::HS_WCAPeriodicCellLists<2> pot_cellA(eps, sca, radii, boxvec, x, rcut, factor);
        pele::HS_WCAPeriodicCellLists<2> pot_cellB(eps, sca, radii, boxvec, x, rcut, factor + 0.2);
        pele::Array<double> g_no_cells(x.size());
        pele::Array<double> g_cellA(x.size());
        pele::Array<double> g_cellB(x.size());
        pele::Array<double> h_no_cells(x.size() * x.size());
        pele::Array<double> h_cellA(h_no_cells.size());
        pele::Array<double> h_cellB(h_no_cells.size());
        const double egh_no_cells = pot_no_cells.get_energy_gradient_hessian(x, g_no_cells, h_no_cells);
        const double egh_cellA = pot_cellA.get_energy_gradient_hessian(x, g_cellA, h_cellA);
        const double egh_cellB = pot_cellB.get_energy_gradient_hessian(x, g_cellB, h_cellB);
        EXPECT_DOUBLE_EQ(e_no_cells, egh_no_cells);
        EXPECT_DOUBLE_EQ(e_no_cells, egh_cellA);
        EXPECT_DOUBLE_EQ(e_no_cells, egh_cellB);
        for (size_t i = 0; i < g_no_cells.size(); ++i) {
            EXPECT_DOUBLE_EQ(g_no_cells[i], g_cellA[i]);
            EXPECT_DOUBLE_EQ(g_no_cells[i], g_cellB[i]);
        }
        for (size_t i = 0; i < h_no_cells.size(); ++i) {
            EXPECT_DOUBLE_EQ(h_no_cells[i], h_cellA[i]);
            EXPECT_DOUBLE_EQ(h_no_cells[i], h_cellB[i]);
        }
    }
}

TEST_F(CellIterTestMoreHS_WCA2D, HSWCAMinimzation_Works) {
    auto pot_no_cells = std::make_shared<pele::HS_WCAPeriodic<2> >(eps, sca, radii, boxvec);
    auto pot_cells = std::make_shared<pele::HS_WCAPeriodicCellLists<2> >(eps, sca, radii, boxvec, x, rcut, 1);
    pele::MODIFIED_FIRE opt_no_cells(pot_no_cells, x, .1, 1, 1);
    pele::MODIFIED_FIRE opt_cells(pot_cells, x, .1, 1, 1);
    opt_no_cells.run();
    opt_cells.run();
    auto x_opt_no_cells = opt_no_cells.get_x();
    auto x_opt_cells = opt_no_cells.get_x();
    const auto e_opt_no_cells = pot_no_cells->get_energy(x_opt_no_cells);
    const auto e_opt_cells = pot_cells->get_energy(x_opt_cells);
    EXPECT_DOUBLE_EQ(e_opt_no_cells, e_opt_cells);
    for (size_t i = 0; i < x_opt_no_cells.size(); ++i) {
        EXPECT_DOUBLE_EQ(x_opt_no_cells[i], x_opt_cells[i]);
    }
}

TEST_F(CellIterTestMoreHS_WCA2D, HSWCAMinimzationPerformance_Changes) {
    auto pot_no_cells = std::make_shared<pele::HS_WCAPeriodic<2> >(eps, sca, radii, boxvec);
    auto pot_cells = std::make_shared<pele::HS_WCAPeriodicCellLists<2> >(eps, sca, radii, boxvec, x, 1, 0.3);
    pele::MODIFIED_FIRE opt_no_cells(pot_no_cells, x, .1, 1, 1);
    pele::MODIFIED_FIRE opt_cells(pot_cells, x, .1, 1, 1);
    const size_t repetitions = 50;
    const auto start_time = clock();
    for (size_t i = 0; i < repetitions; ++i) {
        opt_no_cells.reset(x);
        opt_no_cells.run();
    }
    const auto intermediate_1 = clock();
    for (size_t i = 0; i < repetitions; ++i) {
        opt_cells.reset(x);
        opt_cells.run();
    }
    const auto intermediate_2 = clock();
    auto x_opt_no_cells = opt_no_cells.get_x();
    auto x_opt_cells = opt_no_cells.get_x();
    const auto e_opt_no_cells = pot_no_cells->get_energy(x_opt_no_cells);
    const auto e_opt_cells = pot_cells->get_energy(x_opt_cells);
    EXPECT_DOUBLE_EQ(e_opt_no_cells, e_opt_cells);
    for (size_t i = 0; i < x_opt_no_cells.size(); ++i) {
        EXPECT_DOUBLE_EQ(x_opt_no_cells[i], x_opt_cells[i]);
    }
    const auto time_no_cells = intermediate_1 - start_time;
    const auto time_cells = intermediate_2 - intermediate_1;
    EXPECT_LE(time_cells, time_no_cells);
}
*/
