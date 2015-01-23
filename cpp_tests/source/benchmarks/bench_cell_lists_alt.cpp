#include <random>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "pele/cell_lists.h"
#include "pele/lj_cut.h"
#include "pele/matrix.h"

#include "bench_utils.hpp"

using namespace pele;
using std::string;

struct Rand {
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution;
    Rand()
        : distribution(0,1)
    {}
    double get() { return distribution(generator); }
};

struct LJCellListMaker {
    Rand random_double;
    double rcut;
    double ncellx_scale;

    LJCellListMaker(Rand & rand_, double rcut_, double ncellx_scale_)
        : random_double(rand_),
          rcut(rcut),
          ncellx_scale(ncellx_scale_)
    {}

    std::shared_ptr<BasePotential> get_potential_coords(size_t natoms, Array<double> x)
    {
        double density = 1.2;
        double boxl = std::pow(natoms / density * (4./3 * M_PI), 1./3);
        std::cout << "box length " << boxl << std::endl;
        Array<double> boxvec(3, boxl);

        assert(x.size() == 3*natoms);
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] = random_double.get() * boxl;
        }

        auto lj = std::make_shared<LJCutPeriodicCellLists<3> >(4., 4., rcut, boxvec, ncellx_scale);
        lj->get_energy(x); // this does the initialization

        return lj;
    }
};


double bench_potential(std::shared_ptr<BasePotential> pot, Array<double> x, size_t neval)
{
    auto grad = x.copy();
    Timer t;
    t.start();
    for (size_t i = 0; i < neval; ++i) {
        // change x by some amount and recompute the energy
        double dx = .1;
        if (i % 5 == 0) dx *= -1;
        x[i % x.size()] += dx;
        pot->get_energy_gradient(x, grad);
//        if (i % 500 == 0) {
//            std::cout << i << " energy " << energy << "\n";
//        }
    }
    t.stop();

    return t.get();
}



double benc_potential_natoms(LJCellListMaker & pot_maker,
        size_t natoms,
        size_t neval)
{
    std::cout << "\n";
    std::cout << "======================================================\n";
    std::cout << "benchmarking " << neval << " energy evaluations for " << natoms << " atoms\n";


    Array<double> x(3*natoms);

    Timer timer;
    timer.start();
    auto pot = pot_maker.get_potential_coords(natoms, x);
    timer.stop();
    double time_init = timer.get();
    std::cout <<  "  time to initialize potential: " << time_init << "\n";

    double t = bench_potential(pot, x, neval);
    std::cout <<  "  total time: " << t << "\n";
    std::cout <<  "  timer per evaluation: " << t / neval << "\n";
    std::cout <<  "  timer per evaluation per atom: " << t / neval / natoms << "\n";
    std::cout << "  natoms neval time_init time_eval\n";
    std::cout << "  all times:"
            << " " << natoms
            << " " << neval
            << " " << time_init
            << " " << t
            << "\n";
    return t;
}


int main(int argc, char ** argv)
{
    std::cout << std::setprecision(16);
    Rand r;


    double rcut = 2.;

    double ncellx_scale = 1.;

    double neval = 10000;
    size_t natoms = 10;
    double const target_time_per_run = 5.; // in seconds
    double time_per_eval_per_atom = 1e-5;

    LJCellListMaker pot_maker(r, rcut, ncellx_scale);

    while (neval >= 1) {
        double t = benc_potential_natoms(pot_maker, natoms, std::round(neval));
        time_per_eval_per_atom = t / neval / natoms;
        natoms *= 2;
        neval = target_time_per_run / (time_per_eval_per_atom * natoms);
    }

}
