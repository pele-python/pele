#include "pele/array.h"
#include "pele/vecn.h"
#include "pele/cell_lists.h"
#include "pele/distance.h"
#include "pele/hs_wca.h"
#include "pele/inversepower.h"
#include "pele/lj_cut.h"
#include "pele/modified_fire.h"
#include "test_utils.hpp"

#include <iostream>
#include <stdexcept>
#include <gtest/gtest.h>
#include <random>
#include <ctime>
#include <omp.h>
#include <algorithm>

using pele::Array;
using pele::InversePowerPeriodic;
using pele::InversePower_interaction;

static double const EPS = std::numeric_limits<double>::min();
#define EXPECT_NEAR_RELATIVE(A, B, T)  EXPECT_NEAR(A/(fabs(A)+fabs(B) + EPS), B/(fabs(A)+fabs(B) + EPS), T)

template<size_t ndim>
class stupid_counter {
private:
    std::vector<size_t> m_count;
public:
    stupid_counter()
    {
        #ifdef _OPENMP
        m_count = std::vector<size_t>(omp_get_max_threads(), 0);
        #else
        m_count = std::vector<size_t>(1, 0);
        #endif
    }

    void insert_atom_pair(const size_t atom_i, const size_t atom_j, const size_t isubdom)
    {
        #ifdef _OPENMP
        m_count[omp_get_thread_num()]++;
        #else
        m_count[0]++;
        #endif
    }

    double get_count() {
        return std::accumulate(m_count.begin(), m_count.end(), 0);
    }
 };

template<typename distance_policy>
size_t get_nr_unique_pairs(Array<double> coords, pele::CellLists<distance_policy> & cl)
{
    stupid_counter<distance_policy::_ndim> counter;
    cl.update(coords);
    auto looper = cl.get_atom_pair_looper(counter);
    looper.loop_through_atom_pairs();
    return counter.get_count();
}

template<typename distance_policy>
size_t get_direct_nr_unique_pairs(std::shared_ptr<distance_policy> dist,
        const double max_distance, pele::Array<double> x)
{
    static const size_t m_ndim = distance_policy::_ndim;
    size_t nr_unique_pairs = 0;
    const size_t natoms = x.size() / m_ndim;
    for (size_t i = 0; i < natoms; ++i) {
        for (size_t j = i + 1; j < natoms; ++j) {
            double rij[m_ndim];
            const double* xi = x.data() + i * m_ndim;
            const double* xj = x.data() + j * m_ndim;
            dist->get_rij(rij, xi, xj);
            double r2 = 0;
            for (size_t k = 0; k < m_ndim; ++k) {
                r2 += rij[k] * rij[k];
            }
            nr_unique_pairs += (r2 <= (max_distance * max_distance));
        }
    }
    return nr_unique_pairs;
}


class CellListsTest : public ::testing::Test {
public:
    double pow, eps, etrue, etrue_r, rcut, sca;
    Array<double> x, g, gnum, radii, boxvec;
    void SetUp(){
        #ifdef _OPENMP
        omp_set_num_threads(1);
        #endif
    	pow = 2.5;
    	eps = 1;
    	x = Array<double>(9);
        x[0] = 0.1;
        x[1] = 0.2;
        x[2] = 0.3;
        x[3] = 0.44;
        x[4] = 0.55;
        x[5] = 1.66;
        x[6] = 0.88;
        x[7] = 1.1;
        x[8] = 2.49;
        radii = Array<double>(3);
        boxvec = Array<double>(3, 5);
        for (size_t j = 0; j < 3; ++j) {
            double center = 0;
            for (size_t k = 0; k < 3; ++k) {
                center += x[k * 3 + j] / double(3);
            }
            for (size_t k = 0; k < 3; ++k) {
                x[k * 3 + j] -= center;
            }
        }
        double f = 1.;
        radii[0] = .91 * f;
        radii[1] = 1.1 * f;
        radii[2] = 1.13 * f;
        etrue = 0.03493116137645523;
        g = Array<double>(x.size());
        gnum = Array<double>(x.size());
        sca = 1.2;
        rcut = 2 * (1 + sca) * *std::max_element(radii.begin(), radii.end());
    }
};

//test number of distinguishable pairs
TEST_F(CellListsTest, Number_of_neighbors){
    pele::CellLists<> cell(std::make_shared<pele::periodic_distance<3> >(boxvec), boxvec, boxvec[0]);
    pele::CellLists<> cell2(std::make_shared<pele::periodic_distance<3> >(boxvec), boxvec, boxvec[0], 1);
    pele::CellLists<> cell3(std::make_shared<pele::periodic_distance<3> >(boxvec), boxvec, boxvec[0], 4.2);
    pele::CellLists<> cell4(std::make_shared<pele::periodic_distance<3> >(boxvec), boxvec, boxvec[0], 5);
    size_t count = 3u;
    ASSERT_EQ(count, get_nr_unique_pairs(x, cell));
    ASSERT_EQ(count, get_nr_unique_pairs(x, cell2));
    ASSERT_EQ(count, get_nr_unique_pairs(x, cell3));
    ASSERT_EQ(count, get_nr_unique_pairs(x, cell4));
}

TEST_F(CellListsTest, Number_of_neighbors_Cartesian){
    pele::CellLists<pele::cartesian_distance<3> > cell(std::make_shared<pele::cartesian_distance<3> >(), boxvec, boxvec[0]);
    pele::CellLists<pele::cartesian_distance<3> > cell2(std::make_shared<pele::cartesian_distance<3> >(), boxvec, boxvec[0], 1);
    pele::CellLists<pele::cartesian_distance<3> > cell3(std::make_shared<pele::cartesian_distance<3> >(), boxvec, boxvec[0], 4.2);
    pele::CellLists<pele::cartesian_distance<3> > cell4(std::make_shared<pele::cartesian_distance<3> >(), boxvec, boxvec[0], 5);
    size_t count = 3u;
    ASSERT_EQ(count, get_nr_unique_pairs(x, cell));
    ASSERT_EQ(count, get_nr_unique_pairs(x, cell2));
    ASSERT_EQ(count, get_nr_unique_pairs(x, cell3));
    ASSERT_EQ(count, get_nr_unique_pairs(x, cell4));
}

TEST_F(CellListsTest, NumberNeighborsDifferentRcut_Works){
    auto dist = std::make_shared<pele::periodic_distance<3> >(boxvec);
    pele::CellLists<> cell(dist, boxvec, boxvec[0]);
    pele::CellLists<> cell2(dist, boxvec, boxvec[0], 1);
    pele::CellLists<> cell3(dist, boxvec, boxvec[0], 4.2);
    pele::CellLists<> cell4(dist, boxvec, boxvec[0], 5);
    size_t count_ref = get_direct_nr_unique_pairs(dist, boxvec[0], x);
    size_t count = get_nr_unique_pairs(x, cell);
    size_t count2 = get_nr_unique_pairs(x, cell2);
    size_t count3 = get_nr_unique_pairs(x, cell3);
    size_t count4 = get_nr_unique_pairs(x, cell4);
    ASSERT_EQ(3u, count_ref);
    ASSERT_EQ(count_ref, count);
    ASSERT_EQ(count_ref, count2);
    ASSERT_EQ(count_ref, count3);
    ASSERT_EQ(count_ref, count4);
}

TEST_F(CellListsTest, NumberNeighborsDifferentRcut_WorksLeesEdwards){
    for(double shear = 0.0; shear <= 1.0; shear += 0.01) {
        auto dist = std::make_shared<pele::leesedwards_distance<3> >(boxvec, shear);
        pele::CellLists<pele::leesedwards_distance<3>> cell(dist, boxvec, boxvec[0]);
        pele::CellLists<pele::leesedwards_distance<3>> cell2(dist, boxvec, boxvec[0], 1);
        pele::CellLists<pele::leesedwards_distance<3>> cell3(dist, boxvec, boxvec[0], 4.2);
        pele::CellLists<pele::leesedwards_distance<3>> cell4(dist, boxvec, boxvec[0], 5);
        size_t count_ref = get_direct_nr_unique_pairs(dist, boxvec[0], x);
        size_t count = get_nr_unique_pairs(x, cell);
        size_t count2 = get_nr_unique_pairs(x, cell2);
        size_t count3 = get_nr_unique_pairs(x, cell3);
        size_t count4 = get_nr_unique_pairs(x, cell4);
        ASSERT_EQ(3u, count_ref);
        ASSERT_EQ(count_ref, count);
        ASSERT_EQ(count_ref, count2);
        ASSERT_EQ(count_ref, count3);
        ASSERT_EQ(count_ref, count4);
    }
}

TEST_F(CellListsTest, NumberNeighborsDifferentRcut_WorksCartesian){
    auto dist = std::make_shared<pele::cartesian_distance<3> >();
    pele::CellLists<pele::cartesian_distance<3> > cell(dist, boxvec, boxvec[0]);
    pele::CellLists<pele::cartesian_distance<3> > cell2(dist, boxvec, boxvec[0], 1);
    pele::CellLists<pele::cartesian_distance<3> > cell3(dist, boxvec, boxvec[0], 4.2);
    pele::CellLists<pele::cartesian_distance<3> > cell4(dist, boxvec, boxvec[0], 5);
    size_t count = get_direct_nr_unique_pairs(dist, boxvec[0], x);
    size_t count2 = get_nr_unique_pairs(x, cell2);
    size_t count3 = get_nr_unique_pairs(x, cell3);
    size_t count4 = get_nr_unique_pairs(x, cell4);
    ASSERT_EQ(3u, count);
    ASSERT_EQ(count, count2);
    ASSERT_EQ(count, count3);
    ASSERT_EQ(count, count4);
}

TEST_F(CellListsTest, getEnergy_nan){
    pele::InversePowerPeriodicCellLists<3> pot_cell(pow, eps, radii, boxvec, 1.0);
    x.assign(NAN);
    const double e = pot_cell.get_energy(x);
    ASSERT_TRUE(isnan(e));
}

TEST_F(CellListsTest, getEnergy_inf){
    pele::InversePowerPeriodicCellLists<3> pot_cell(pow, eps, radii, boxvec, 1.0);
    x.assign(INFINITY);
    const double e = pot_cell.get_energy(x);
    ASSERT_TRUE(isnan(e));
}

TEST_F(CellListsTest, getEnergyGradient_nan){
    pele::InversePowerPeriodicCellLists<3> pot_cell(pow, eps, radii, boxvec, 1.0);
    x.assign(NAN);
    pele::Array<double> g(x.size());
    g.assign(NAN);
    const double e = pot_cell.get_energy_gradient(x, g);
    ASSERT_TRUE(isnan(e));
    ASSERT_TRUE(std::all_of(g.begin(), g.end(),
                        [](double elem) { return isnan(elem); }
                            ));
}

TEST_F(CellListsTest, getEnergyGradient_inf){
    pele::InversePowerPeriodicCellLists<3> pot_cell(pow, eps, radii, boxvec, 1.0);
    x.assign(INFINITY);
    pele::Array<double> g(x.size());
    g.assign(INFINITY);
    const double e = pot_cell.get_energy_gradient(x, g);
    ASSERT_TRUE(isnan(e));
    ASSERT_TRUE(std::all_of(g.begin(), g.end(),
                        [](double elem) { return isnan(elem); }
                            ));
}

TEST_F(CellListsTest, getEnergyGradientHessian_nan){
    pele::InversePowerPeriodicCellLists<3> pot_cell(pow, eps, radii, boxvec, 1.0);
    x.assign(NAN);
    pele::Array<double> g(x.size());
    g.assign(NAN);
    pele::Array<double> h(x.size() * x.size());
    h.assign(NAN);
    const double e = pot_cell.get_energy_gradient_hessian(x, g, h);
    ASSERT_TRUE(isnan(e));
    ASSERT_TRUE(std::all_of(g.begin(), g.end(),
                        [](double elem) { return isnan(elem); }
                            ));
    ASSERT_TRUE(std::all_of(h.begin(), h.end(),
                        [](double elem) { return isnan(elem); }
                            ));
}

TEST_F(CellListsTest, getEnergyGradientHessian_inf){
    pele::InversePowerPeriodicCellLists<3> pot_cell(pow, eps, radii, boxvec, 1.0);
    x.assign(INFINITY);
    pele::Array<double> g(x.size());
    g.assign(INFINITY);
    pele::Array<double> h(x.size() * x.size());
    h.assign(INFINITY);
    const double e = pot_cell.get_energy_gradient_hessian(x, g, h);
    ASSERT_TRUE(isnan(e));
    ASSERT_TRUE(std::all_of(g.begin(), g.end(),
                        [](double elem) { return isnan(elem); }
                            ));
    ASSERT_TRUE(std::all_of(h.begin(), h.end(),
                        [](double elem) { return isnan(elem); }
                            ));
}

TEST_F(CellListsTest, Energy_Works){
    pele::InversePowerPeriodicCellLists<3> pot_cell(pow, eps, radii, boxvec, 1.0);
    pele::InversePowerPeriodicCellLists<3> pot_cell2(pow, eps, radii, boxvec, 2.0);
    pele::InversePowerPeriodicCellLists<3> pot_cell3(pow, eps, radii, boxvec, 3.0);
    pele::InversePowerPeriodicCellLists<3> pot_cell4(pow, eps, radii, boxvec, 4.0);
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

TEST_F(CellListsTest, ChangeCoords_Works){
    pele::InversePowerPeriodicCellLists<3> pot_cell(pow, eps, radii, boxvec, .1);
    pele::InversePowerPeriodic<3> pot(pow, eps, radii, boxvec);
    double ecell = pot_cell.get_energy(x);
    double etrue = pot.get_energy(x);
    ASSERT_NEAR(ecell, etrue, 1e-10);

    for (size_t i = 0; i < x.size(); ++i) {
        x[i] += i;
    }
    ecell = pot_cell.get_energy(x);
    etrue = pot.get_energy(x);
    ASSERT_NEAR(ecell, etrue, 1e-10);

    for (size_t i = 0; i < x.size(); ++i) {
        x[i] += (i+4)*4;
    }
    ecell = pot_cell.get_energy(x);
    etrue = pot.get_energy(x);
    ASSERT_NEAR(ecell, etrue, 1e-10);
}


//TEST_F(CellListsTest, EnergyCartesian_Works){
//    pele::InversePowerCellLists<3> pot_cell(pow, eps, radii, boxvec, rcut, 1.0);
//    pele::InversePowerCellLists<3> pot_cell2(pow, eps, radii, boxvec, rcut, 2.0);
//    pele::InversePowerCellLists<3> pot_cell3(pow, eps, radii, boxvec, rcut, 3.0);
//    pele::InversePowerCellLists<3> pot_cell4(pow, eps, radii, boxvec, rcut, 4.0);
//    pele::InversePower<3> pot(pow, eps, radii);
//    const double ecell = pot_cell.get_energy(x);
//    const double ecell2 = pot_cell2.get_energy(x);
//    const double ecell3 = pot_cell3.get_energy(x);
//    const double ecell4 = pot_cell4.get_energy(x);
//    const double etrue = pot.get_energy(x);
//    ASSERT_NEAR(ecell, etrue, 1e-10);
//    ASSERT_NEAR(ecell2, etrue, 1e-10);
//    ASSERT_NEAR(ecell3, etrue, 1e-10);
//    ASSERT_NEAR(ecell4, etrue, 1e-10);
//}

TEST_F(CellListsTest, EnergyGradient_AgreesWithNumerical){
    pele::InversePowerPeriodic<3> pot_no_cells(pow, eps, radii, boxvec);
    const double etrue = pot_no_cells.get_energy(x);
    const size_t N = 3;
    std::vector<std::shared_ptr<pele::InversePowerPeriodicCellLists<3> > > pot;
    for (size_t i = 0; i < N; ++i) {
        pot.push_back(std::make_shared<pele::InversePowerPeriodicCellLists<3> >(
                pow, eps, radii, boxvec, 1 + i));
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

//TEST_F(CellListsTest, EnergyGradientCartesian_AgreesWithNumerical){
//    pele::InversePower<3> pot_no_cells(pow, eps, radii);
//    const double etrue = pot_no_cells.get_energy(x);
//    const size_t N = 3;
//    std::vector<std::shared_ptr<pele::InversePowerCellLists<3> > > pot;
//    for (size_t i = 0; i < N; ++i) {
//        pot.push_back(std::make_shared<pele::InversePowerCellLists<3> >(
//                pow, eps, radii, boxvec, rcut, 1 + i));
//    }
//    pot.swap(pot);
//    std::vector<double> e(N, 0);
//    std::vector<double> ecomp(N, 0);
//    for (size_t i = 0; i < N; ++i) {
//        e.at(i) = pot.at(i)->get_energy_gradient(x, g);
//        ecomp.at(i) = pot.at(i)->get_energy(x);
//        pot.at(i)->numerical_gradient(x, gnum, 1e-6);
//        for (size_t k = 0; k < 6; ++k) {
//            ASSERT_NEAR(g[k], gnum[k], 1e-6);
//        }
//    }
//    for (size_t i = 0; i < N; ++i) {
//        ASSERT_NEAR(e.at(i), ecomp.at(i), 1e-10);
//        ASSERT_NEAR(e.at(i), etrue, 1e-10);
//    }
//}

TEST_F(CellListsTest, EnergyGradientHessian_AgreesWithNumerical){
    pele::InversePowerPeriodic<3> pot_no_cells(pow, eps, radii, boxvec);
    const double etrue = pot_no_cells.get_energy(x);
    Array<double> g_no_cells(x.size()) ;
    Array<double> h_no_cells(x.size() * x.size());
    pot_no_cells.get_energy_gradient_hessian(x, g_no_cells, h_no_cells);
    for (size_t i = 0; i < 3; ++i) {
        pele::InversePowerPeriodicCellLists<3> pot(pow, eps, radii, boxvec, 1.0 + i);
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

//TEST_F(CellListsTest, EnergyGradientHessianCartesian_AgreesWithNumerical){
//    pele::InversePower<3> pot_no_cells(pow, eps, radii);
//    const double etrue = pot_no_cells.get_energy(x);
//    Array<double> g_no_cells(x.size()) ;
//    Array<double> h_no_cells(x.size() * x.size());
//    pot_no_cells.get_energy_gradient_hessian(x, g_no_cells, h_no_cells);
//    for (size_t i = 0; i < 3; ++i) {
//        pele::InversePowerCellLists<3> pot(pow, eps, radii, boxvec, rcut, 1.0 + i);
//        Array<double> h(x.size() * x.size());
//        Array<double> hnum(h.size());
//        const double e = pot.get_energy_gradient_hessian(x, g, h);
//        const double ecomp = pot.get_energy(x);
//        pot.numerical_gradient(x, gnum);
//        pot.numerical_hessian(x, hnum);
//        EXPECT_NEAR(e, ecomp, 1e-10);
//        EXPECT_NEAR(etrue, ecomp, 1e-10);
//        for (size_t i = 0; i < g.size(); ++i) {
//            ASSERT_NEAR(g[i], gnum[i], 1e-10);
//            ASSERT_NEAR(g[i], g_no_cells[i], 1e-10);
//        }
//        for (size_t i = 0; i < h.size(); ++i) {
//            ASSERT_NEAR(h[i], hnum[i], 1e-10);
//            ASSERT_NEAR(h[i], h_no_cells[i], 1e-10);
//        }
//    }
//}

TEST_F(CellListsTest, HS_WCAEnergy_Works){
    pele::HS_WCAPeriodicCellLists<3> pot_cell(eps, sca, radii, boxvec, 1);
    pele::HS_WCAPeriodicCellLists<3> pot_cell2(eps, sca, radii, boxvec, 1.1);
    pele::HS_WCAPeriodicCellLists<3> pot_cell3(eps, sca, radii, boxvec, 1.2);
    pele::HS_WCAPeriodicCellLists<3> pot_cell4(eps, sca, radii, boxvec, 1.3);
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

TEST_F(CellListsTest, HS_WCAEnergyLeesEdwards_Works){
    for(double shear = 0.0; shear <= 1.0; shear += 0.01) {
        pele::HS_WCALeesEdwardsCellLists<3> pot_cell(eps, sca, radii, boxvec, shear, 1);
        pele::HS_WCALeesEdwardsCellLists<3> pot_cell2(eps, sca, radii, boxvec, shear, 1.1);
        pele::HS_WCALeesEdwardsCellLists<3> pot_cell3(eps, sca, radii, boxvec, shear, 1.2);
        pele::HS_WCALeesEdwardsCellLists<3> pot_cell4(eps, sca, radii, boxvec, shear, 1.3);
        pele::HS_WCALeesEdwards<3> pot_no_cells(eps, sca, radii, boxvec, shear);
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
}

TEST_F(CellListsTest, HS_WCAEnergyCartesian_Works){
    pele::HS_WCACellLists<3> pot_cell(eps, sca, radii, boxvec, 1);
    pele::HS_WCACellLists<3> pot_cell2(eps, sca, radii, boxvec, 1.1);
    pele::HS_WCACellLists<3> pot_cell3(eps, sca, radii, boxvec, 1.2);
    pele::HS_WCACellLists<3> pot_cell4(eps, sca, radii, boxvec, 1.3);
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

TEST_F(CellListsTest, Norm_SimpleTest){
    pele::HS_WCACellLists<3, 6> pot(eps, sca, radii, boxvec, 1);

    pele::Array<double> x = {0.0, 1.0, 1.0, 3.0, 4.0, 0.0, 2.0, 2.0, 2.0, 3.0, 0.0, 2.0, 4.0, 1.0, 2.0, 3.0, 3.0, 2.0};
    double norm = pot.compute_norm(x);
    EXPECT_DOUBLE_EQ(norm, 5.0);
}


class CellListsTestHomogeneous3D : public ::testing::Test {
public:
    size_t nparticles;
    size_t boxdim;
    double boxlength;
    pele::Array<double> boxvec;
    size_t ndof;
    size_t seed;
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> distribution;
    pele::Array<double> x;
    void SetUp(){
        #ifdef _OPENMP
        omp_set_num_threads(1);
        #endif
        nparticles = 200;
        boxdim = 3;
        boxlength = 42;
        boxvec = pele::Array<double>(boxdim, boxlength);
        ndof = nparticles * boxdim;
        seed = 42;
        generator = std::mt19937_64(seed);
        distribution = std::uniform_real_distribution<double>(-0.5 * boxlength, 0.5 * boxlength);
        x = pele::Array<double>(ndof);
        for (size_t i = 0; i < ndof; ++i) {
            x[i] = distribution(generator);
        }
    }
};

TEST_F(CellListsTestHomogeneous3D, GridAndSpacing_Works) {
    pele::CellLists<> cell_one(std::make_shared<pele::periodic_distance<3> >(boxvec), boxvec, boxvec[0]);
    cell_one.update(x);
    EXPECT_EQ(cell_one.get_nr_cells(), 1u);
    EXPECT_EQ(cell_one.get_nr_cellsx(), 1u);
    //std::cout << "nr_unique_pairs: one:\n" << get_nr_unique_pairs(x, cell_one) << "\n";
    pele::CellLists<> cell_two(std::make_shared<pele::periodic_distance<3> >(boxvec), boxvec, boxvec[0] / 2);
    cell_two.update(x);
    EXPECT_EQ(cell_two.get_nr_cells(), 8u);
    EXPECT_EQ(cell_two.get_nr_cellsx(), 2u);
    //std::cout << "nr_unique_pairs: two:\n" << get_nr_unique_pairs(x, cell_two) << "\n";
    pele::CellLists<> cell_three(std::make_shared<pele::periodic_distance<3> >(boxvec), boxvec, boxvec[0] / 3);
    cell_three.update(x);
    EXPECT_EQ(cell_three.get_nr_cells(), 27u);
    EXPECT_EQ(cell_three.get_nr_cellsx(), 3u);
}

TEST_F(CellListsTestHomogeneous3D, GridAndSpacingLeesEdwards_Works) {
    for(double shear = 0.0; shear <= 1.0; shear += 0.01) {
        pele::CellLists<pele::leesedwards_distance<3> > cell_one(
            std::make_shared<pele::leesedwards_distance<3> >(boxvec, shear), boxvec, boxvec[0]);
        cell_one.update(x);
        EXPECT_EQ(cell_one.get_nr_cells(), 1u);
        EXPECT_EQ(cell_one.get_nr_cellsx(), 1u);
        //std::cout << "nr_unique_pairs: one:\n" << get_nr_unique_pairs(x, cell_one) << "\n";
        pele::CellLists<pele::leesedwards_distance<3> > cell_two(
            std::make_shared<pele::leesedwards_distance<3> >(boxvec, shear), boxvec, boxvec[0] / 2);
        cell_two.update(x);
        EXPECT_EQ(cell_two.get_nr_cells(), 8u);
        EXPECT_EQ(cell_two.get_nr_cellsx(), 2u);
        //std::cout << "nr_unique_pairs: two:\n" << get_nr_unique_pairs(x, cell_two) << "\n";
        pele::CellLists<pele::leesedwards_distance<3> > cell_three(
            std::make_shared<pele::leesedwards_distance<3> >(boxvec, shear), boxvec, boxvec[0] / 3);
        cell_three.update(x);
        EXPECT_EQ(cell_three.get_nr_cells(), 27u);
        EXPECT_EQ(cell_three.get_nr_cellsx(), 3u);
    }
}

TEST_F(CellListsTestHomogeneous3D, GridAndSpacingCartesian_Works) {
    pele::CellLists<pele::cartesian_distance<3> > cell_one(std::make_shared<pele::cartesian_distance<3> >(), boxvec, boxvec[0]);
    EXPECT_EQ(cell_one.get_nr_cells(), 1u);
    EXPECT_EQ(cell_one.get_nr_cellsx(), 1u);
    //std::cout << "nr_unique_pairs: one:\n" << get_nr_unique_pairs(x, cell_one) << "\n";
    pele::CellLists<pele::cartesian_distance<3> > cell_two(std::make_shared<pele::cartesian_distance<3> >(), boxvec, boxvec[0] / 2);
    EXPECT_EQ(cell_two.get_nr_cells(), 8u);
    EXPECT_EQ(cell_two.get_nr_cellsx(), 2u);
    //std::cout << "nr_unique_pairs: two:\n" << get_nr_unique_pairs(x, cell_two) << "\n";
    pele::CellLists<pele::cartesian_distance<3> > cell_three(std::make_shared<pele::cartesian_distance<3> >(), boxvec, boxvec[0] / 3);
    EXPECT_EQ(cell_three.get_nr_cells(), 27u);
    EXPECT_EQ(cell_three.get_nr_cellsx(), 3u);
}


class CellListsTestHomogeneous2D : public ::testing::Test {
public:
    size_t nparticles;
    size_t boxdim;
    double boxlength;
    pele::Array<double> boxvec;
    size_t ndof;
    size_t seed;
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> distribution;
    pele::Array<double> x;
    void SetUp(){
        #ifdef _OPENMP
        omp_set_num_threads(1);
        #endif
        nparticles = 200;
        boxdim = 2;
        boxlength = 42;
        boxvec = pele::Array<double>(boxdim, boxlength);
        ndof = nparticles * boxdim;
        seed = 42;
        generator = std::mt19937_64(seed);
        distribution = std::uniform_real_distribution<double>(-0.5 * boxlength, 0.5 * boxlength);
        x = pele::Array<double>(ndof);
        for (size_t i = 0; i < ndof; ++i) {
            x[i] = distribution(generator);
        }
    }
};

TEST_F(CellListsTestHomogeneous2D, GridAndSpacing_Works) {
    pele::CellLists<pele::periodic_distance<2> > cell_one(std::make_shared<pele::periodic_distance<2> >(boxvec), boxvec, boxvec[0]);
    EXPECT_EQ(cell_one.get_nr_cells(), 1u);
    EXPECT_EQ(cell_one.get_nr_cellsx(), 1u);
    //std::cout << "nr_unique_pairs: one:\n" << get_nr_unique_pairs(x, cell_one) << "\n";
    pele::CellLists<pele::periodic_distance<2> > cell_two(std::make_shared<pele::periodic_distance<2> >(boxvec), boxvec, boxvec[0] / 2);
    EXPECT_EQ(cell_two.get_nr_cells(), 4u);
    EXPECT_EQ(cell_two.get_nr_cellsx(), 2u);
    //std::cout << "nr_unique_pairs: two:\n" << get_nr_unique_pairs(x, cell_two) << "\n";
    pele::CellLists<pele::periodic_distance<2> > cell_three(std::make_shared<pele::periodic_distance<2> >(boxvec), boxvec, boxvec[0] / 3);
    EXPECT_EQ(cell_three.get_nr_cells(), 9u);
    EXPECT_EQ(cell_three.get_nr_cellsx(), 3u);
}

TEST_F(CellListsTestHomogeneous2D, GridAndSpacingCartesian_Works) {
    pele::CellLists<pele::cartesian_distance<2> > cell_one(std::make_shared<pele::cartesian_distance<2> >(), boxvec, boxvec[0]);
    EXPECT_EQ(cell_one.get_nr_cells(), 1u);
    EXPECT_EQ(cell_one.get_nr_cellsx(), 1u);
    //std::cout << "nr_unique_pairs: one:\n" << get_nr_unique_pairs(x, cell_one) << "\n";
    pele::CellLists<pele::cartesian_distance<2> > cell_two(std::make_shared<pele::cartesian_distance<2> >(), boxvec, boxvec[0] / 2);
    EXPECT_EQ(cell_two.get_nr_cells(), 4u);
    EXPECT_EQ(cell_two.get_nr_cellsx(), 2u);
    //std::cout << "nr_unique_pairs: two:\n" << get_nr_unique_pairs(x, cell_two) << "\n";
    pele::CellLists<pele::cartesian_distance<2> > cell_three(std::make_shared<pele::cartesian_distance<2> >(), boxvec, boxvec[0] / 3);
    EXPECT_EQ(cell_three.get_nr_cells(), 9u);
    EXPECT_EQ(cell_three.get_nr_cellsx(), 3u);
}


class LatticeNeighborsTest : public ::testing::Test {
public:
    static const size_t ndim = 3;
    typedef pele::periodic_distance<ndim> dist_t;
    Array<double> boxvec;
    double rcut;
    Array<size_t> ncells_vec;
    std::shared_ptr<dist_t> dist;

    virtual void SetUp(){
        #ifdef _OPENMP
        omp_set_num_threads(1);
        #endif
        boxvec = Array<double>(3, 10);
        boxvec[1] += 1;
        boxvec[2] += 2;
        rcut = 20.; // large rcut means all cells are neighbors
        dist = std::make_shared<dist_t> (boxvec);
        ncells_vec = Array<size_t>(ndim);
        ncells_vec[0] = 2;
        ncells_vec[1] = 4;
        ncells_vec[2] = 20;
    }
};

TEST_F(LatticeNeighborsTest, LargeRcut_Works)
{
    static size_t const ndim = 3;
    typedef pele::periodic_distance<ndim> dist_t;

    pele::LatticeNeighbors<dist_t> lattice(dist, boxvec, rcut, ncells_vec);

    size_t icell = 8+3;
    auto v = lattice.global_ind_to_cell_vec(icell);
    // std::cout << v << std::endl;
    ASSERT_EQ(icell, lattice.cell_vec_to_global_ind(v));

    icell = 47;
    v = lattice.global_ind_to_cell_vec(icell);
    // std::cout << v << std::endl;
    ASSERT_EQ(icell, lattice.cell_vec_to_global_ind(v));


    auto neibs = lattice.find_all_global_neighbor_inds(0);
    ASSERT_EQ(neibs.size(), lattice.m_ncells);

    // check there are no duplicates
    std::set<size_t> s(neibs.begin(), neibs.end());
    ASSERT_EQ(neibs.size(), s.size());

    std::vector< std::vector< std::array<long*, 2> > > pairs_inner(lattice.m_nsubdoms);
    std::vector< std::vector< std::array<long*, 2> > > pairs_boundary(lattice.m_nsubdoms);
    std::vector< std::vector<long> > cells(lattice.m_nsubdoms);
    for (size_t isubdom = 0; isubdom < lattice.m_nsubdoms; isubdom++) {
        cells[isubdom] = std::vector<long>(lattice.cell_vec_to_global_ind(ncells_vec) / lattice.m_nsubdoms);
    }
    lattice.find_neighbor_pairs(pairs_inner, pairs_boundary, cells);
    size_t count_neighbors = 0;
    for (size_t isubdom = 0; isubdom < lattice.m_nsubdoms; isubdom++) {
        count_neighbors += pairs_inner[isubdom].size() + pairs_boundary[isubdom].size();
    }
    ASSERT_EQ(count_neighbors, lattice.m_ncells * (lattice.m_ncells+1)/2);

    pele::VecN<ndim> xpos(0.1);
    size_t icell1, isubdom1;
    lattice.position_to_local_ind(xpos.data(), icell1, isubdom1);

    xpos[0] += 2.;
    xpos[1] += 3.4;
    xpos[2] += .9;
    size_t icell2, isubdom2;
    lattice.position_to_local_ind(xpos.data(), icell2, isubdom2);
    ASSERT_NE(icell1, icell2);
}

TEST_F(LatticeNeighborsTest, SmallRcut_Works2)
{
    static size_t const ndim = 3;
    typedef pele::periodic_distance<ndim> dist_t;
    rcut = .1; // small rcut means only adjacent cells are neighbors

    pele::LatticeNeighbors<dist_t> lattice(dist, boxvec, rcut, ncells_vec);

    auto neibs = lattice.find_all_global_neighbor_inds(0);
    ASSERT_EQ(neibs.size(), size_t(2*3*3));

    // check there are no duplicates
    std::set<size_t> s(neibs.begin(), neibs.end());
    ASSERT_EQ(neibs.size(), s.size());

    std::vector< std::vector< std::array<long*, 2> > > pairs_inner(lattice.m_nsubdoms);
    std::vector< std::vector< std::array<long*, 2> > > pairs_boundary(lattice.m_nsubdoms);
    std::vector< std::vector<long> > cells(lattice.m_nsubdoms);
    for (size_t isubdom = 0; isubdom < lattice.m_nsubdoms; isubdom++) {
        cells[isubdom] = std::vector<long>(lattice.cell_vec_to_global_ind(ncells_vec) / lattice.m_nsubdoms);
    }
    lattice.find_neighbor_pairs(pairs_inner, pairs_boundary, cells);
    size_t count_neighbors = 0;
    for (size_t isubdom = 0; isubdom < lattice.m_nsubdoms; isubdom++) {
        count_neighbors += pairs_inner[isubdom].size() + pairs_boundary[isubdom].size();
    }
    ASSERT_EQ(count_neighbors, (neibs.size() - 1) * lattice.m_ncells / 2 + lattice.m_ncells );
}

TEST_F(LatticeNeighborsTest, NonPeriodic_Works2)
{
    static size_t const ndim = 3;
    double rcut = .1; // small rcut means only adjacent cells are neighbors
    typedef pele::cartesian_distance<ndim> dist_t;
    auto cart_dist = std::make_shared<dist_t> ();

    pele::LatticeNeighbors<dist_t> lattice(cart_dist, boxvec, rcut, ncells_vec);

    auto neibs = lattice.find_all_global_neighbor_inds(0);
    ASSERT_EQ(neibs.size(), size_t(2*2*2));

    // check there are no duplicates
    std::set<size_t> s(neibs.begin(), neibs.end());
    ASSERT_EQ(neibs.size(), s.size());

    neibs = lattice.find_all_global_neighbor_inds(2);
    ASSERT_EQ(neibs.size(), size_t(2*3*2));
}

TEST_F(LatticeNeighborsTest, positionToCellVec_BoxBoundaryWorks)
{
    static size_t const ndim = 3;
    typedef pele::periodic_distance<ndim> dist_t;
    pele::LatticeNeighbors<dist_t> lattice(dist, boxvec, rcut, ncells_vec);

    double boxboundary = boxvec[0] * 0.5;
    for (int i = -20; i <= 20; i++) {
        for (int j = -20; j <= 20; j++) {
            pele::Array<double> coords(ndim, 0);
            coords[0] = i * boxboundary + j * std::numeric_limits<double>::epsilon();
            pele::VecN<ndim, size_t> cell_vec = lattice.position_to_cell_vec(coords.data());
            EXPECT_LT(cell_vec[0], lattice.m_ncells_vec[0]);
        }
    }
}

using pele::SimplePairwisePotential;
using pele::lj_interaction_cut_smooth;
using pele::periodic_distance;

template<size_t NDIM>
class LJCutPeriodicN : public SimplePairwisePotential<lj_interaction_cut_smooth, periodic_distance<NDIM> > {
public:
    LJCutPeriodicN(double C6, double C12, double rcut, Array<double> const boxvec)
        : SimplePairwisePotential< lj_interaction_cut_smooth, periodic_distance<NDIM> > (
                std::make_shared<lj_interaction_cut_smooth>(C6, C12, rcut),
                std::make_shared<periodic_distance<NDIM> >(boxvec)
                )
    {}
};

template<size_t ndim>
void test_rectangular_cell_lists(const double rcut, Array<double> boxvec, Array<double> x)
{
    periodic_distance<ndim> dist(boxvec);
    std::shared_ptr<pele::BasePotential> pot_c = std::make_shared<pele::LJCutPeriodicCellLists<ndim> >(4., 4., rcut, boxvec, 1);
    std::shared_ptr<pele::BasePotential> pot_p = std::make_shared<LJCutPeriodicN<ndim> >(4., 4., rcut, boxvec);
    EXPECT_DOUBLE_EQ(pot_c->get_energy(x), pot_p->get_energy(x));
    Array<double> g_c(x.size());
    auto g_p = g_c.copy();
    Array<double> h_c(x.size() * x.size());
    auto h_p = h_c.copy();
    pot_p->get_energy_gradient(x, g_p);
    pot_c->get_energy_gradient(x, g_c);
    for (size_t i = 0; i < x.size(); ++i) {
        EXPECT_NEAR_RELATIVE(g_p[i], g_c[i], 1e-10);
    }
    pot_p->get_energy_gradient_hessian(x, g_p, h_p);
    pot_c->get_energy_gradient_hessian(x, g_c, h_c);
    for (size_t i = 0; i < g_p.size(); ++i) {
        EXPECT_NEAR_RELATIVE(g_p[i], g_c[i], 1e-10);
    }
    for (size_t i = 0; i < h_p.size(); ++i) {
        EXPECT_NEAR_RELATIVE(h_p[i], h_c[i], 1e-10);
    }
    std::shared_ptr<pele::GradientOptimizer> optimizer_c = std::make_shared<pele::MODIFIED_FIRE>(pot_c, x, 1, 1, .1);
    std::shared_ptr<pele::GradientOptimizer> optimizer_p = std::make_shared<pele::MODIFIED_FIRE>(pot_p, x, 1, 1, .1);
    optimizer_c->run();
    optimizer_p->run();
    auto xo_p = optimizer_p->get_x();
    auto xo_c = optimizer_c->get_x();
    EXPECT_NEAR(pot_c->get_energy(xo_c), pot_p->get_energy(xo_p), 1e-10);
    for (size_t i = 0; i < xo_p.size(); ++i) {
        EXPECT_NEAR_RELATIVE(xo_p[i], xo_c[i], 1e-9);
    }
}

void get_boxvec(const double rcut, const size_t ndim, Array<double>& boxvec)
{
    boxvec = Array<double>(ndim);
    for (size_t i = 0; i < ndim; ++i) {
        boxvec[i] = std::pow((2 * rcut + 1), i + 1) + 0.8 * i;
    }
}

template<size_t ndim>
void get_boxvec_x0L(const double rcut, const size_t N, Array<double>& boxvec, Array<double>& x)
{
    get_boxvec(rcut, ndim, boxvec);
    x = Array<double>(ndim * N);
    std::mt19937_64 gen(42);
    for (size_t i = 0; i < N; ++i) {
        for (size_t k = 0; k < ndim; ++k) {
            x[i * ndim + k] = std::uniform_real_distribution<double>{0, boxvec[k]}(gen);
        }
    }
}

template<size_t ndim>
void get_boxvec_x05(const double rcut, const size_t N, Array<double>& boxvec, Array<double>& x)
{
    get_boxvec(rcut, ndim, boxvec);
    x = Array<double>(ndim * N);
    std::mt19937_64 gen(44);
    for (size_t i = 0; i < N; ++i) {
        for (size_t k = 0; k < ndim; ++k) {
            x[i * ndim + k] = std::uniform_real_distribution<double>{-0.5 * boxvec[k], 0.5 * boxvec[k]}(gen);
        }
    }
}

template<size_t ndim>
void get_boxvec_x_all(const double rcut, const size_t N, Array<double>& boxvec, Array<double>& x)
{
    get_boxvec(rcut, ndim, boxvec);
    x = Array<double>(ndim * N);
    std::mt19937_64 gen(46);
    for (size_t i = 0; i < ndim * N; ++i) {
        x[i] = std::uniform_real_distribution<double>{-1e10, 1e10}(gen);
    }
}

TEST_F(LatticeNeighborsTest, RectangleWorks)
{
    rcut = 2.5;
    const size_t N = 25;
    Array<double> x;
    get_boxvec_x0L<2>(rcut, N, boxvec, x);
    test_rectangular_cell_lists<2>(rcut, boxvec, x);
    get_boxvec_x0L<3>(rcut, N, boxvec, x);
    test_rectangular_cell_lists<3>(rcut, boxvec, x);
    get_boxvec_x05<2>(rcut, N, boxvec, x);
    test_rectangular_cell_lists<2>(rcut, boxvec, x);
    get_boxvec_x05<3>(rcut, N, boxvec, x);
    test_rectangular_cell_lists<3>(rcut, boxvec, x);
    get_boxvec_x_all<2>(rcut, N, boxvec, x);
    test_rectangular_cell_lists<2>(rcut, boxvec, x);
}


class OpenMPCellListsTest : public ::testing::Test {
public:
    size_t seed;
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> distribution;
    size_t nparticles;
    size_t ndim;
    size_t ndof;
    double eps;
    double sca;
    double r_hs;
    Array<double> x;
    Array<double> radii;
    Array<double> boxvec;
    double rcut;
    virtual void SetUp(){
        #ifndef _OPENMP
        throw std::runtime_error("This test only works with OpenMP enabled compilation");
        #endif
        seed = 42;
        generator = std::mt19937_64(seed);
        distribution = std::uniform_real_distribution<double>(-1, 1);
        nparticles = 50;
        ndim = 2;
        ndof = nparticles * ndim;
        eps = 1;
        r_hs = 1;
        x = Array<double>(ndof);
        radii = Array<double>(nparticles);
        for (size_t i = 0; i < nparticles; ++i) {
            radii[i] = r_hs;
        }
        sca = 0.2;
        rcut = 2 * (1 + sca) * *std::max_element(radii.begin(), radii.end());
    }

    void create_coords() {
        // Order atoms like a stair
        Array<double> coords(2);
        size_t k = 0;
        for(int i = 0; i < nparticles; i++) {
            coords[k] += (1 + (0.6 + 0.5 * distribution(generator)) * sca) * 2 * r_hs;
            x[2*i] = coords[0];
            x[2*i + 1] = coords[1];
            if(k == 0) {
                k = 1;
            } else {
                k = 0;
            }
        }
        boxvec = Array<double>(ndim, std::max<double>(fabs(*std::max_element(x.data(), x.data() + ndof)), fabs(*std::min_element(x.data(), x.data() + ndof))) + rcut);
    }
};

TEST_F(OpenMPCellListsTest, HSWCAEnergyLeesEdwards_Works) {
    for(size_t nthreads = 2; nthreads <= 4; nthreads++) {
        #ifdef _OPENMP
        omp_set_num_threads(nthreads);
        #endif
        for(double shear = 0.0; shear <= 0.2; shear += 0.1) {
            create_coords();
            pele::HS_WCALeesEdwards<2> pot_no_cells(eps, sca, radii, boxvec, shear);
            pele::HS_WCALeesEdwardsCellLists<2> pot_cell(eps, sca, radii, boxvec, shear, 1.0);
            for(size_t new_coords = 0; new_coords < 10; new_coords++) {
                create_coords();
                const double e_no_cells = pot_no_cells.get_energy(x);
                for(size_t rep_same = 0; rep_same < 100; rep_same++) {
                    const double e_cell = pot_cell.get_energy(x);
                    EXPECT_TRUE(almostEqual(e_no_cells, e_cell, 8));
                }
            }
        }
    }
}

TEST_F(OpenMPCellListsTest, HSWCAEnergyGradientLeesEdwards_Works) {
    for(size_t nthreads = 2; nthreads <= 4; nthreads++) {
        #ifdef _OPENMP
        omp_set_num_threads(nthreads);
        #endif
        for(double shear = 0.0; shear <= 0.2; shear += 0.1) {
            create_coords();

            pele::HS_WCALeesEdwards<2> pot_no_cells(eps, sca, radii, boxvec, shear);
            pele::HS_WCALeesEdwardsCellLists<2> pot_cell(eps, sca, radii, boxvec, shear, 1.0);
            for(size_t new_coords = 0; new_coords < 10; new_coords++) {
                create_coords();

                const double e_no_cells = pot_no_cells.get_energy(x);
                pele::Array<double> g_no_cells(x.size());
                pele::Array<double> g_cell(x.size());
                const double eg_no_cells = pot_no_cells.get_energy_gradient(x, g_no_cells);
                for(size_t rep_same = 0; rep_same < 100; rep_same++) {
                    const double eg_cell = pot_cell.get_energy_gradient(x, g_cell);
                    EXPECT_TRUE(almostEqual(e_no_cells, eg_no_cells, 8));
                    EXPECT_TRUE(almostEqual(e_no_cells, eg_cell, 8));
                    for (size_t i = 0; i < g_no_cells.size(); ++i) {
                        EXPECT_TRUE(almostEqual(g_no_cells[i], g_cell[i], 8));
                    }
                }
            }
        }
    }
}

TEST_F(OpenMPCellListsTest, HSWCAEnergyGradientHessianLeesEdwards_Works) {
    for(size_t nthreads = 2; nthreads <= 4; nthreads++) {
        #ifdef _OPENMP
        omp_set_num_threads(nthreads);
        #endif
        for(double shear = 0.0; shear <= 0.2; shear += 0.1) {
            create_coords();
            pele::HS_WCALeesEdwards<2> pot_no_cells(eps, sca, radii, boxvec, shear);
            pele::HS_WCALeesEdwardsCellLists<2> pot_cell(eps, sca, radii, boxvec, shear, 1.0);
            for(size_t new_coords = 0; new_coords < 10; new_coords++) {
                create_coords();
                const double e_no_cells = pot_no_cells.get_energy(x);
                pele::Array<double> g_no_cells(x.size());
                pele::Array<double> g_cell(x.size());
                pele::Array<double> h_no_cells(x.size() * x.size());
                pele::Array<double> h_cell(h_no_cells.size());
                const double egh_no_cells = pot_no_cells.get_energy_gradient_hessian(x, g_no_cells, h_no_cells);
                for(size_t rep_same = 0; rep_same < 100; rep_same++) {
                    const double egh_cell = pot_cell.get_energy_gradient_hessian(x, g_cell, h_cell);
                    EXPECT_TRUE(almostEqual(e_no_cells, egh_no_cells, 8));
                    EXPECT_TRUE(almostEqual(e_no_cells, egh_cell, 8));
                    for (size_t i = 0; i < g_no_cells.size(); ++i) {
                        EXPECT_TRUE(almostEqual(g_no_cells[i], g_cell[i], 8));
                    }
                    for (size_t i = 0; i < h_no_cells.size(); ++i) {
                        EXPECT_TRUE(almostEqual(h_no_cells[i], h_cell[i], 8));
                    }
                }
            }
        }
    }
}
