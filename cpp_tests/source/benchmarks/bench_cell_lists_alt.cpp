#include <random>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <ctime>

#include "pele/cell_lists.h"
#include "pele/lj_cut.h"
#include "pele/matrix.h"

#include "bench_utils.hpp"

using namespace pele;
using std::string;

class Timer {
public:
    double tstart, tstop;

    void start() { tstart = clock(); }
    void stop() { tstop = clock(); }
    double get() { return (tstop - tstart) / CLOCKS_PER_SEC; }
};

int main(int argc, char ** argv)
{
    std::cout << std::setprecision(16);
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0,1);

    size_t natoms = 1600;
    double rcut = 2.;
    double density = 1.2;
    double boxl = std::pow(natoms / density * (4./3 * M_PI), 1./3);
    std::cout << "box length " << boxl << std::endl;
    Array<double> boxvec(3, boxl);

    Array<double> x(3*natoms);
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = distribution(generator) * boxl;
    }

    if (argc > 1) {
        auto fname = string(argv[1]);
        std::cout << "reading from file " << fname << std::endl;
        x = coords_from_file(fname);
    }

    // print the coords array array
//    pele::MatrixAdapter<double> m(x, 3);
//    std::cout << "initial coordinates\n";
//    std::cout << m << std::endl;

    double ncellx_scale = 1.;
    auto lj = std::make_shared<LJCutPeriodicCellLists<3> >(4., 4., rcut, boxvec, ncellx_scale);

    std::cout << "energy " << lj->get_energy(x) << std::endl;

    Timer timer;
    auto grad = x.copy();
    size_t neval = 10000;
    timer.start();
    for (size_t i = 0; i < neval; ++i) {
        // change x by some amount and recompute the energy
        double dx = .1;
        if (i % 5 == 0) dx *= -1;
        x[i % x.size()] += dx;
        double energy = lj->get_energy_gradient(x, grad);
        if (i % 500 == 0) {
            std::cout << i << " energy " << energy << "\n";
        }
    }
    timer.stop();
    std::cout <<  "total time: " << timer.get() << "\n";
    std::cout <<  "timer per evaluation: " << timer.get() / neval << "\n";
    std::cout <<  "timer per evaluation per atom: " << timer.get() / neval / natoms << "\n";


}
