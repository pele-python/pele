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


struct LJCellListMaker {
    MyRNG & random_double;
    double rcut;
    double ncellx_scale;

    LJCellListMaker(MyRNG & rand_, double rcut_, double ncellx_scale_)
        : random_double(rand_),
          rcut(rcut_),
          ncellx_scale(ncellx_scale_)
    {}

    std::shared_ptr<BasePotential> get_potential_coords(Array<double> x)
    {
        double density = 1.2;
        size_t natoms = x.size() / 3;
        double boxl = std::pow(natoms / density * (4./3 * M_PI), 1./3);
        std::cout << "box length " << boxl << std::endl;
        Array<double> boxvec(3, boxl);

        assert(x.size() == 3*natoms);
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] = random_double.get() * boxl;
        }

        std::cout << rcut << " " << boxvec << " " << ncellx_scale << "\n";
        auto lj = std::make_shared<LJCutPeriodicCellLists<3> >(4., 4., rcut, boxvec, ncellx_scale);
        double energy = lj->get_energy(x); // this does the initialization
        std::cout << "initial energy " << energy << "\n";
        return lj;
    }
};

int main(int argc, char ** argv)
{
    std::cout << std::setprecision(16);
    MyRNG r;


    double rcut = 2.;

    double ncellx_scale = 1.;

    double neval = 10000;
    size_t natoms = 10;
    double const target_time_per_run = 5.; // in seconds

    LJCellListMaker pot_maker(r, rcut, ncellx_scale);

    while (neval >= 1) {
        double t = bench_potential_natoms(pot_maker, natoms, std::round(neval));
        double time_per_eval_per_atom = t / neval / natoms;
        natoms *= 2;
        neval = target_time_per_run / (time_per_eval_per_atom * natoms);
    }

}
