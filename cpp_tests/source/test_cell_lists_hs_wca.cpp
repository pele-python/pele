#include <omp.h>
#include <gtest/gtest.h>

#include "pele/array.h"
#include "pele/cell_lists.h"
#include "pele/distance.h"
#include "pele/hs_wca.h"
#include "pele/inversepower.h"
#include "pele/modified_fire.h"
#include "test_utils.hpp"

using pele::Array;

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

class CellListsTestMoreHS_WCA : public ::testing::Test {
public:
    double pow;
    size_t seed;
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> distribution;
    size_t nparticles;
    size_t ndim;
    size_t ndof;
    double eps;
    double sca;
    Array<double> x;
    Array<double> g;
    Array<double> gnum;
    Array<double> radii;
    Array<double> boxvec;
    double rcut;
    void SetUp(){
        #ifdef _OPENMP
        omp_set_num_threads(1);
        #endif
        pow = 2.5;
        seed = 42;
        generator = std::mt19937_64(seed);
        distribution = std::uniform_real_distribution<double>(0, 0.05);
        nparticles = 5;
        ndim = 3;
        ndof = nparticles * ndim;
        eps = 1;
        x = Array<double>(ndof);
        for (size_t k = 0; k < ndof; ++k) {
            x[k] = double(k + 1) / double(10) + 0.2 * (k / nparticles) * (1 + distribution(generator));
        }
        for (size_t j = 0; j < ndim; ++j) {
            double center = 0;
            for (size_t k = 0; k < nparticles; ++k) {
                center += x[k * ndim + j] / static_cast<double>(nparticles);
            }
            for (size_t k = 0; k < nparticles; ++k) {
                x[k * ndim + j] -= center;
            }
        }
        radii = Array<double>(nparticles);
        for (size_t i = 0; i < nparticles; ++i) {
            radii[i] = (0.1 + distribution(generator));
        }
        g = Array<double>(x.size());
        gnum = Array<double>(x.size());
        sca = 1.2;
        rcut = 2 * (1 + sca) * *std::max_element(radii.data(), radii.data() + nparticles);
        //std::cout << "rcut: " << rcut << std::endl;
        const double L = 2 * std::max<double>(fabs(*std::max_element(x.data(), x.data() + ndof)), fabs(*std::min_element(x.data(), x.data() + ndof))) + 2 * rcut;
        //std::cout << "L: " << L << std::endl;
        boxvec = Array<double>(ndim, L);
    }
};

TEST_F(CellListsTestMoreHS_WCA, Number_of_neighbors){
    pele::CellLists<> cell(std::make_shared<pele::periodic_distance<3> >(boxvec), boxvec, boxvec[0]);
    pele::CellLists<> cell2(std::make_shared<pele::periodic_distance<3> >(boxvec), boxvec, boxvec[0], 1);
    pele::CellLists<> cell3(std::make_shared<pele::periodic_distance<3> >(boxvec), boxvec, boxvec[0], 2);
    pele::CellLists<> cell4(std::make_shared<pele::periodic_distance<3> >(boxvec), boxvec, boxvec[0], 4);
    size_t count = nparticles * (nparticles - 1) / 2;
    ASSERT_EQ(count, get_nr_unique_pairs(x, cell));
    ASSERT_EQ(count, get_nr_unique_pairs(x, cell2));
    ASSERT_EQ(count, get_nr_unique_pairs(x, cell3));
    ASSERT_EQ(count, get_nr_unique_pairs(x, cell4));
}

TEST_F(CellListsTestMoreHS_WCA, Number_of_neighbors_LeesEdwards)
{
    for(double shear = 0.0; shear <= 1.0; shear += 0.01) {
        pele::CellLists<pele::leesedwards_distance<3> > cell(std::make_shared<pele::leesedwards_distance<3> >(boxvec, shear), boxvec, boxvec[0]);
        pele::CellLists<pele::leesedwards_distance<3> > cell2(std::make_shared<pele::leesedwards_distance<3> >(boxvec, shear), boxvec, boxvec[0], 1);
        pele::CellLists<pele::leesedwards_distance<3> > cell3(std::make_shared<pele::leesedwards_distance<3> >(boxvec, shear), boxvec, boxvec[0], 2);
        pele::CellLists<pele::leesedwards_distance<3> > cell4(std::make_shared<pele::leesedwards_distance<3> >(boxvec, shear), boxvec, boxvec[0], 4);
        size_t count = nparticles * (nparticles - 1) / 2;
        ASSERT_EQ(count, get_nr_unique_pairs(x, cell));
        ASSERT_EQ(count, get_nr_unique_pairs(x, cell2));
        ASSERT_EQ(count, get_nr_unique_pairs(x, cell3));
        ASSERT_EQ(count, get_nr_unique_pairs(x, cell4));
    }
}

TEST_F(CellListsTestMoreHS_WCA, Number_of_neighbors_Cartesian){
    for (size_t ii = 0; ii < ndof; ++ii) {
        EXPECT_LE(x[ii], 0.5 * boxvec[0]);
        EXPECT_LE(-0.5 * boxvec[0], x[ii]);
    }
    for (size_t i = 0; i < radii.size(); ++i) {
        EXPECT_LE(0, radii[i]);
    }
    pele::CellLists<pele::cartesian_distance<3> > cell(std::make_shared<pele::cartesian_distance<3> >(), boxvec, boxvec[0]);
    size_t count = get_nr_unique_pairs(x, cell);
    ASSERT_EQ(nparticles * (nparticles - 1) / 2, count);
//    ASSERT_EQ(count, static_cast<unsigned int>(cell.end() - cell.begin()));
    ASSERT_EQ(count, get_nr_unique_pairs(x, cell));
}


TEST_F(CellListsTestMoreHS_WCA, EnergyMoreParticles_Works){
    for (size_t ii = 0; ii < ndof; ++ii) {
        EXPECT_LE(x[ii], 0.5 * boxvec[0]);
        EXPECT_LE(-0.5 * boxvec[0], x[ii]);
    }
    pele::InversePowerPeriodicCellLists<3> pot_cell(pow, eps, radii, boxvec, .1);
    pele::InversePowerPeriodicCellLists<3> pot_cell2(pow, eps, radii, boxvec, .2);
    pele::InversePowerPeriodicCellLists<3> pot_cell3(pow, eps, radii, boxvec, .3);
    pele::InversePowerPeriodicCellLists<3> pot_cell_(pow, eps, radii, boxvec, .1);
    pele::InversePowerPeriodicCellLists<3> pot_cell2_(pow, eps, radii, boxvec, .2);
    pele::InversePowerPeriodicCellLists<3> pot_cell3_(pow, eps, radii, boxvec, .3);
    pele::InversePowerPeriodic<3> pot(pow, eps, radii, boxvec);
    const double ecell = pot_cell.get_energy(x);
    const double ecell2 = pot_cell2.get_energy(x);
    const double ecell3 = pot_cell3.get_energy(x);
    const double etrue = pot.get_energy(x);
    const double ecell_ = pot_cell_.get_energy(x);
    const double ecell2_ = pot_cell2_.get_energy(x);
    const double ecell3_ = pot_cell3_.get_energy(x);
    EXPECT_DOUBLE_EQ(ecell, etrue);
    EXPECT_DOUBLE_EQ(ecell2, etrue);
    EXPECT_DOUBLE_EQ(ecell3, etrue);
    EXPECT_DOUBLE_EQ(ecell_, etrue);
    EXPECT_DOUBLE_EQ(ecell2_, etrue);
    EXPECT_DOUBLE_EQ(ecell3_, etrue);
}

//TEST_F(CellListsTestMoreHS_WCA, EnergyMoreParticlesCartesian_Works){
//    for (size_t ii = 0; ii < ndof; ++ii) {
//        EXPECT_LE(x[ii], 0.5 * boxvec[0]);
//        EXPECT_LE(-0.5 * boxvec[0], x[ii]);
//    }
//    pele::InversePowerCellLists<3> pot_cell(pow, eps, radii, boxvec, rcut, .1);
//    pele::InversePowerCellLists<3> pot_cell2(pow, eps, radii, boxvec, rcut, .2);
//    pele::InversePowerCellLists<3> pot_cell3(pow, eps, radii, boxvec, rcut, .3);
//    pele::InversePowerCellLists<3> pot_cell_(pow, eps, radii, boxvec, boxvec[0], 1);
//    pele::InversePowerCellLists<3> pot_cell2_(pow, eps, radii, boxvec, boxvec[0], 2);
//    pele::InversePowerCellLists<3> pot_cell3_(pow, eps, radii, boxvec, boxvec[0], 3);
//    pele::InversePower<3> pot(pow, eps, radii);
//    const double ecell = pot_cell.get_energy(x);
//    const double ecell2 = pot_cell2.get_energy(x);
//    const double ecell3 = pot_cell3.get_energy(x);
//    const double etrue = pot.get_energy(x);
//    const double ecell_ = pot_cell_.get_energy(x);
//    const double ecell2_ = pot_cell2_.get_energy(x);
//    const double ecell3_ = pot_cell3_.get_energy(x);
//    EXPECT_DOUBLE_EQ(ecell, etrue);
//    EXPECT_DOUBLE_EQ(ecell2, etrue);
//    EXPECT_DOUBLE_EQ(ecell3, etrue);
//    EXPECT_DOUBLE_EQ(ecell_, etrue);
//    EXPECT_DOUBLE_EQ(ecell2_, etrue);
//    EXPECT_DOUBLE_EQ(ecell3_, etrue);
//}

TEST_F(CellListsTestMoreHS_WCA, HSWCAEnergy_Works) {
    for (size_t ii = 0; ii < ndof; ++ii) {
        EXPECT_LE(x[ii], 0.5 * boxvec[0]);
        EXPECT_LE(-0.5 * boxvec[0], x[ii]);
    }
    pele::HS_WCAPeriodic<3> pot_no_cells(eps, sca, radii, boxvec);
    const double e_no_cells = pot_no_cells.get_energy(x);
    for (size_t factor = 1; factor < 2; ++factor) {
        pele::HS_WCAPeriodicCellLists<3> pot_cellA(eps, sca, radii, boxvec, factor * 0.01);
        pele::HS_WCAPeriodicCellLists<3> pot_cellB(eps, sca, radii, boxvec, (factor + 0.01) * 0.01);
        const double e_cellA = pot_cellA.get_energy(x);
        const double e_cellB = pot_cellB.get_energy(x);
        EXPECT_DOUBLE_EQ(e_no_cells, e_cellA);
        EXPECT_DOUBLE_EQ(e_no_cells, e_cellB);
    }
}

TEST_F(CellListsTestMoreHS_WCA, HSWCAEnergyLeesEdwards_Works) {
    for(double shear = 0.0; shear <= 1.0; shear += 0.01) {
        for (size_t ii = 0; ii < ndof; ++ii) {
            EXPECT_LE(x[ii], 0.5 * boxvec[0]);
            EXPECT_LE(-0.5 * boxvec[0], x[ii]);
        }
        pele::HS_WCA<3> pot_no_cells(eps, sca, radii);
        const double e_no_cells = pot_no_cells.get_energy(x);
        for (size_t factor = 1; factor < 2; ++factor) {
            //std::cout << "factor: " << factor << std::endl;
            pele::HS_WCALeesEdwardsCellLists<3> pot_cell1(eps, sca, radii, boxvec, shear);
            const double e_cell1 = pot_cell1.get_energy(x);
            pele::HS_WCALeesEdwardsCellLists<3> pot_cellA(eps, sca, radii, boxvec, shear, factor);
            pele::HS_WCALeesEdwardsCellLists<3> pot_cellB(eps, sca, radii, boxvec, shear);
            const double e_cellA = pot_cellA.get_energy(x);
            const double e_cellB = pot_cellB.get_energy(x);
            //std::cout << "e_no_cells: " << e_no_cells << "\n";
            //std::cout << "e_cell1: " << e_cell1 << std::endl;
            EXPECT_DOUBLE_EQ(e_cell1, e_cellA);
            EXPECT_DOUBLE_EQ(e_cellA, e_cellB);
            EXPECT_DOUBLE_EQ(e_no_cells, e_cellA);
            EXPECT_DOUBLE_EQ(e_no_cells, e_cellB);
        }
    }
}

TEST_F(CellListsTestMoreHS_WCA, HSWCAEnergyCartesian_Works) {
    for (size_t ii = 0; ii < ndof; ++ii) {
        EXPECT_LE(x[ii], 0.5 * boxvec[0]);
        EXPECT_LE(-0.5 * boxvec[0], x[ii]);
    }
    pele::HS_WCA<3> pot_no_cells(eps, sca, radii);
    const double e_no_cells = pot_no_cells.get_energy(x);
    for (size_t factor = 1; factor < 2; ++factor) {
        //std::cout << "factor: " << factor << std::endl;
        pele::HS_WCACellLists<3> pot_cell1(eps, sca, radii, boxvec);
        const double e_cell1 = pot_cell1.get_energy(x);
        pele::HS_WCACellLists<3> pot_cellA(eps, sca, radii, boxvec, factor);
        pele::HS_WCACellLists<3> pot_cellB(eps, sca, radii, boxvec);
        const double e_cellA = pot_cellA.get_energy(x);
        const double e_cellB = pot_cellB.get_energy(x);
        //std::cout << "e_no_cells: " << e_no_cells << "\n";
        //std::cout << "e_cell1: " << e_cell1 << std::endl;
        EXPECT_DOUBLE_EQ(e_cell1, e_cellA);
        EXPECT_DOUBLE_EQ(e_cellA, e_cellB);
        EXPECT_DOUBLE_EQ(e_no_cells, e_cellA);
        EXPECT_DOUBLE_EQ(e_no_cells, e_cellB);
    }
}

TEST_F(CellListsTestMoreHS_WCA, HSWCAEnergyGradient_Works) {
    for (size_t ii = 0; ii < ndof; ++ii) {
        EXPECT_LE(x[ii], 0.5 * boxvec[0]);
        EXPECT_LE(-0.5 * boxvec[0], x[ii]);
    }
    pele::HS_WCAPeriodic<3> pot_no_cells(eps, sca, radii, boxvec);
    const double e_no_cells = pot_no_cells.get_energy(x);
    for (size_t factor = 1; factor < 3; ++factor) {
        pele::HS_WCAPeriodicCellLists<3> pot_cellA(eps, sca, radii, boxvec, (factor) * 0.2);
        pele::HS_WCAPeriodicCellLists<3> pot_cellB(eps, sca, radii, boxvec, (factor + 0.2) * 0.2);
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

TEST_F(CellListsTestMoreHS_WCA, HSWCAEnergyGradientLeesEdwards_Works) {
    for(double shear = 0.0; shear <= 1.0; shear += 0.01) {
        for (size_t ii = 0; ii < ndof; ++ii) {
            EXPECT_LE(x[ii], 0.5 * boxvec[0]);
            EXPECT_LE(-0.5 * boxvec[0], x[ii]);
        }
        pele::HS_WCALeesEdwards<3> pot_no_cells(eps, sca, radii, boxvec, shear);
        const double e_no_cells = pot_no_cells.get_energy(x);
        for (size_t factor = 1; factor < 3; ++factor) {
            pele::HS_WCALeesEdwardsCellLists<3> pot_cellA(eps, sca, radii, boxvec, shear, (factor) * 0.2);
            pele::HS_WCALeesEdwardsCellLists<3> pot_cellB(eps, sca, radii, boxvec, shear, (factor + 0.2) * 0.2);
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
}

TEST_F(CellListsTestMoreHS_WCA, HSWCAEnergyGradientCartesian_Works) {
    for (size_t ii = 0; ii < ndof; ++ii) {
        EXPECT_LE(x[ii], 0.5 * boxvec[0]);
        EXPECT_LE(-0.5 * boxvec[0], x[ii]);
    }
    pele::HS_WCA<3> pot_no_cells(eps, sca, radii);
    const double e_no_cells = pot_no_cells.get_energy(x);
    for (size_t factor = 1; factor < 3; ++factor) {
        pele::HS_WCACellLists<3> pot_cellA(eps, sca, radii, boxvec, (factor) * 0.2);
        pele::HS_WCACellLists<3> pot_cellB(eps, sca, radii, boxvec, (factor + 0.2) * 0.2);
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

TEST_F(CellListsTestMoreHS_WCA, HSWCAEnergyGradientHessian_Works) {
    for (size_t ii = 0; ii < ndof; ++ii) {
        EXPECT_LE(x[ii], 0.5 * boxvec[0]);
        EXPECT_LE(-0.5 * boxvec[0], x[ii]);
    }
    pele::HS_WCAPeriodic<3> pot_no_cells(eps, sca, radii, boxvec);
    const double e_no_cells = pot_no_cells.get_energy(x);
    for (size_t factor = 1; factor < 3; ++factor) {
        pele::HS_WCAPeriodicCellLists<3> pot_cellA(eps, sca, radii, boxvec, (factor) * 0.2);
        pele::HS_WCAPeriodicCellLists<3> pot_cellB(eps, sca, radii, boxvec, (factor + 0.2) * 0.2);
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

TEST_F(CellListsTestMoreHS_WCA, HSWCAEnergyGradientHessianLeesEdwards_Works) {
    for(double shear = 0.0; shear <= 1.0; shear += 0.01) {
        for (size_t ii = 0; ii < ndof; ++ii) {
            EXPECT_LE(x[ii], 0.5 * boxvec[0]);
            EXPECT_LE(-0.5 * boxvec[0], x[ii]);
        }
        pele::HS_WCALeesEdwards<3> pot_no_cells(eps, sca, radii, boxvec, shear);
        const double e_no_cells = pot_no_cells.get_energy(x);
        for (size_t factor = 1; factor < 3; ++factor) {
            pele::HS_WCALeesEdwardsCellLists<3> pot_cellA(eps, sca, radii, boxvec, shear, (factor) * 0.2);
            pele::HS_WCALeesEdwardsCellLists<3> pot_cellB(eps, sca, radii, boxvec, shear, (factor + 0.2) * 0.2);
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
}

TEST_F(CellListsTestMoreHS_WCA, HSWCAEnergyGradientHessianCartesian_Works) {
    for (size_t ii = 0; ii < ndof; ++ii) {
        EXPECT_LE(x[ii], 0.5 * boxvec[0]);
        EXPECT_LE(-0.5 * boxvec[0], x[ii]);
    }
    pele::HS_WCA<3> pot_no_cells(eps, sca, radii);
    const double e_no_cells = pot_no_cells.get_energy(x);
    for (size_t factor = 1; factor < 3; ++factor) {
        pele::HS_WCACellLists<3> pot_cellA(eps, sca, radii, boxvec, (factor) * 0.2);
        pele::HS_WCACellLists<3> pot_cellB(eps, sca, radii, boxvec, (factor + 0.2) * 0.2);
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

TEST_F(CellListsTestMoreHS_WCA, HSWCAMinimization_Works) {
    for (size_t ii = 0; ii < ndof; ++ii) {
        EXPECT_LE(x[ii], 0.5 * boxvec[0]);
        EXPECT_LE(-0.5 * boxvec[0], x[ii]);
    }
    auto pot_no_cells = std::make_shared<pele::HS_WCAPeriodic<3> >(eps, sca, radii, boxvec);
    auto pot_cells = std::make_shared<pele::HS_WCAPeriodicCellLists<3> >(eps, sca, radii, boxvec, 0.2);
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

TEST_F(CellListsTestMoreHS_WCA, HSWCAMinimizationLeesEdwards_Works) {
    for(double shear = 0.0; shear <= 1.0; shear += 0.01) {
        for (size_t ii = 0; ii < ndof; ++ii) {
            EXPECT_LE(x[ii], 0.5 * boxvec[0]);
            EXPECT_LE(-0.5 * boxvec[0], x[ii]);
        }
        auto pot_no_cells = std::make_shared<pele::HS_WCALeesEdwards<3> >(eps, sca, radii, boxvec, shear);
        auto pot_cells = std::make_shared<pele::HS_WCALeesEdwardsCellLists<3> >(eps, sca, radii, boxvec, shear, 0.2);
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
}

TEST_F(CellListsTestMoreHS_WCA, HSWCAMinimizationCartesian_Works) {
    for (size_t ii = 0; ii < ndof; ++ii) {
        EXPECT_LE(x[ii], 0.5 * boxvec[0]);
        EXPECT_LE(-0.5 * boxvec[0], x[ii]);
    }
    auto pot_no_cells = std::make_shared<pele::HS_WCA<3> >(eps, sca, radii);
    auto pot_cells = std::make_shared<pele::HS_WCACellLists<3> >(eps, sca, radii, boxvec, 0.2);
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

class CellListsTestMoreHS_WCA2D : public ::testing::Test {
public:
    size_t seed;
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> distribution;
    size_t nparticles;
    size_t ndim;
    size_t ndof;
    double eps;
    double sca;
    Array<double> x;
    Array<double> g;
    Array<double> gnum;
    Array<double> radii;
    Array<double> boxvec;
    double rcut;
    virtual void SetUp(){
        #ifdef _OPENMP
        omp_set_num_threads(1);
        #endif
        seed = 42;
        generator = std::mt19937_64(seed);
        distribution = std::uniform_real_distribution<double>(0, 0.05);
        nparticles = 5;
        ndim = 2;
        ndof = nparticles * ndim;
        eps = 1;
        x = Array<double>(ndof);
        for (size_t k = 0; k < ndof; ++k) {
            x[k] = double(k + 1) / double(10) + 0.3 * (k / nparticles) * (1 + distribution(generator));
            x[k] -= distribution(generator) * 0.5;
        }
        radii = Array<double>(nparticles);
        for (size_t i = 0; i < nparticles; ++i) {
            radii[i] = (0.04 + distribution(generator));
        }
        g = Array<double>(x.size());
        gnum = Array<double>(x.size());
        sca = 1.2;
        for (size_t j = 0; j < ndim; ++j) {
            double center = 0;
            for (size_t k = 0; k < nparticles; ++k) {
                center += x[k * ndim + j] / static_cast<double>(nparticles);
            }
            for (size_t k = 0; k < nparticles; ++k) {
                x[k * ndim + j] -= center;
            }
        }
        rcut = 2 * (1 + sca) * *std::max_element(radii.data(), radii.data() + nparticles);
        boxvec = Array<double>(ndim, 2 * std::max<double>(fabs(*std::max_element(x.data(), x.data() + ndof)), fabs(*std::min_element(x.data(), x.data() + ndof))) + rcut);
    }
};

TEST_F(CellListsTestMoreHS_WCA2D, Number_of_neighbors){
    pele::CellLists<pele::periodic_distance<2> > cell(std::make_shared<pele::periodic_distance<2> >(boxvec), boxvec, boxvec[0]);
    pele::CellLists<pele::periodic_distance<2> > cell2(std::make_shared<pele::periodic_distance<2> >(boxvec), boxvec, boxvec[0], 1);
    pele::CellLists<pele::periodic_distance<2> > cell3(std::make_shared<pele::periodic_distance<2> >(boxvec), boxvec, boxvec[0], 3);
    pele::CellLists<pele::periodic_distance<2> > cell4(std::make_shared<pele::periodic_distance<2> >(boxvec), boxvec, boxvec[0], 5);
    size_t count = (nparticles * (nparticles - 1) / 2);
    ASSERT_EQ(count, get_nr_unique_pairs(x, cell));
    ASSERT_EQ(count, get_nr_unique_pairs(x, cell2));
    ASSERT_EQ(count, get_nr_unique_pairs(x, cell3));
    ASSERT_EQ(count, get_nr_unique_pairs(x, cell4));
}

TEST_F(CellListsTestMoreHS_WCA2D, HSWCAEnergy_Works) {
    pele::HS_WCAPeriodic<2> pot_no_cells(eps, sca, radii, boxvec);
    const double e_no_cells = pot_no_cells.get_energy(x);
    for (size_t factor = 1; factor < 3; ++factor) {
        pele::HS_WCAPeriodicCellLists<2> pot_cellA(eps, sca, radii, boxvec, factor);
        pele::HS_WCAPeriodicCellLists<2> pot_cellB(eps, sca, radii, boxvec, factor + 0.2);
        const double e_cellA = pot_cellA.get_energy(x);
        const double e_cellB = pot_cellB.get_energy(x);
        EXPECT_DOUBLE_EQ(e_no_cells, e_cellA);
        EXPECT_DOUBLE_EQ(e_no_cells, e_cellB);
    }
}

TEST_F(CellListsTestMoreHS_WCA2D, HSWCAEnergyLeesEdwards_Works) {
    for(double shear = 0.0; shear <= 1.0; shear += 0.01) {
        pele::HS_WCALeesEdwards<2> pot_no_cells(eps, sca, radii, boxvec, shear);
        const double e_no_cells = pot_no_cells.get_energy(x);
        for (size_t factor = 1; factor < 3; ++factor) {
            pele::HS_WCALeesEdwardsCellLists<2> pot_cellA(eps, sca, radii, boxvec, shear, factor);
            pele::HS_WCALeesEdwardsCellLists<2> pot_cellB(eps, sca, radii, boxvec, shear, factor + 0.2);
            const double e_cellA = pot_cellA.get_energy(x);
            const double e_cellB = pot_cellB.get_energy(x);
            EXPECT_DOUBLE_EQ(e_no_cells, e_cellA);
            EXPECT_DOUBLE_EQ(e_no_cells, e_cellB);
        }
    }
}

TEST_F(CellListsTestMoreHS_WCA2D, HSWCAEnergyCartesian_Works) {
    pele::HS_WCA<2> pot_no_cells(eps, sca, radii);
    pele::HS_WCAPeriodic<2> pot_no_cells_periodic(eps, sca, radii, boxvec);
    const double e_no_cells = pot_no_cells.get_energy(x);
    //const double e_no_cells_periodic = pot_no_cells_periodic.get_energy(x);
    for (size_t factor = 1; factor < 3; ++factor) {
        pele::HS_WCACellLists<2> pot_cellA(eps, sca, radii, boxvec, factor);
        pele::HS_WCAPeriodicCellLists<2> pot_cellA_per(eps, sca, radii, boxvec, factor);
        pele::HS_WCACellLists<2> pot_cellB(eps, sca, radii, boxvec, factor + 0.2);
        const double e_cellA = pot_cellA.get_energy(x);
        const double e_cellB = pot_cellB.get_energy(x);
        //std::cout << "factor: " << factor << "\n";
        //std::cout << "e_no_cells_periodic: " << e_no_cells_periodic << std::endl;
        //std::cout << "radii" << std::endl;
        /*
        for (size_t i = 0; i < radii.size(); ++i) {
            std::cout << radii[i] << "\n";
        }
        std::cout << "boxvec" << std::endl;
        for (size_t i = 0; i < boxvec.size(); ++i) {
            std::cout << boxvec[i] << "\n";
        }
        std::cout << "coords" << std::endl;
        for (size_t i = 0; i < nparticles; ++i) {
            std::cout << x[i * 2] << "\t" << x[i * 2 + 1] << "\t" << radii[i] << "\t" << radii[i] * (1 + sca) << "\n";
        }
        std::cout << "coords_end" << std::endl;
        */
        EXPECT_DOUBLE_EQ(e_no_cells, e_cellA);
        EXPECT_DOUBLE_EQ(e_no_cells, e_cellB);
    }
}

TEST_F(CellListsTestMoreHS_WCA2D, HSWCAEnergyGradient_Works) {
    pele::HS_WCAPeriodic<2> pot_no_cells(eps, sca, radii, boxvec);
    const double e_no_cells = pot_no_cells.get_energy(x);
    for (size_t factor = 1; factor < 3; ++factor) {
        pele::HS_WCAPeriodicCellLists<2> pot_cellA(eps, sca, radii, boxvec, factor);
        pele::HS_WCAPeriodicCellLists<2> pot_cellB(eps, sca, radii, boxvec, factor + 0.2);
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

TEST_F(CellListsTestMoreHS_WCA2D, HSWCAEnergyGradientCartesian_Works) {
    pele::HS_WCA<2> pot_no_cells(eps, sca, radii);
    const double e_no_cells = pot_no_cells.get_energy(x);
    for (size_t factor = 1; factor < 3; ++factor) {
        pele::HS_WCACellLists<2> pot_cellA(eps, sca, radii, boxvec, factor);
        pele::HS_WCACellLists<2> pot_cellB(eps, sca, radii, boxvec, factor + 0.2);
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

TEST_F(CellListsTestMoreHS_WCA2D, HSWCAEnergyGradientHessian_Works) {
    pele::HS_WCAPeriodic<2> pot_no_cells(eps, sca, radii, boxvec);
    const double e_no_cells = pot_no_cells.get_energy(x);
    for (size_t factor = 1; factor < 3; ++factor) {
        pele::HS_WCAPeriodicCellLists<2> pot_cellA(eps, sca, radii, boxvec, factor);
        pele::HS_WCAPeriodicCellLists<2> pot_cellB(eps, sca, radii, boxvec, factor + 0.2);
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

TEST_F(CellListsTestMoreHS_WCA2D, HSWCAEnergyGradientHessianLeesEdwards_Works) {
    for(double shear = 0.0; shear <= 1.0; shear += 0.01) {
        pele::HS_WCALeesEdwards<2> pot_no_cells(eps, sca, radii, boxvec, shear);
        const double e_no_cells = pot_no_cells.get_energy(x);
        for (size_t factor = 1; factor < 3; ++factor) {
            pele::HS_WCALeesEdwardsCellLists<2> pot_cellA(eps, sca, radii, boxvec, shear, factor);
            pele::HS_WCALeesEdwardsCellLists<2> pot_cellB(eps, sca, radii, boxvec, shear, factor + 0.2);
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
}

TEST_F(CellListsTestMoreHS_WCA2D, HSWCAEnergyGradientHessianCartesian_Works) {
    pele::HS_WCA<2> pot_no_cells(eps, sca, radii);
    const double e_no_cells = pot_no_cells.get_energy(x);
    for (size_t factor = 1; factor < 3; ++factor) {
        pele::HS_WCACellLists<2> pot_cellA(eps, sca, radii, boxvec, factor);
        pele::HS_WCACellLists<2> pot_cellB(eps, sca, radii, boxvec, factor + 0.2);
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

TEST_F(CellListsTestMoreHS_WCA2D, HSWCAMinimzation_Works) {
    auto pot_no_cells = std::make_shared<pele::HS_WCAPeriodic<2> >(eps, sca, radii, boxvec);
    auto pot_cells = std::make_shared<pele::HS_WCAPeriodicCellLists<2> >(eps, sca, radii, boxvec, 1);
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
