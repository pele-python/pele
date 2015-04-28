#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <ctime>
#include <random>


#include "pele/array.h"
#include "pele/lj_cut.h"
#include "pele/matrix.h"

class Timer {
public:
    double tstart, tstop;

    void start() { tstart = clock(); }
    void stop() { tstop = clock(); }
    double get() { return (tstop - tstart) / CLOCKS_PER_SEC; }
};


std::vector<double> vector_from_file(std::string fname)
{
    std::vector<double> x;
    std::ifstream fin;
    fin.open(fname.c_str());
    // Prepare a pair of iterators to read the data from cin
    std::istream_iterator<double> eos;
    std::istream_iterator<double> input(fin);
    // No loop is necessary, because you can use copy()
    std::copy(input, eos, std::back_inserter(x));

    fin.close();
    std::cout << "size " << x.size() << std::endl;
    return x;
}


pele::Array<double> coords_from_file(std::string fname)
{
    auto xvec = vector_from_file(fname);
    pele::Array<double> x(xvec);
    return x.copy();
}

struct MyRNG {
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution;
    MyRNG()
        : distribution(0,1)
    {}
    double get() { return distribution(generator); }
};

double bench_potential(std::shared_ptr<pele::BasePotential> pot, pele::Array<double> x, size_t neval)
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


template <class PotentialMaker>
double bench_potential_natoms(PotentialMaker & pot_maker,
        size_t natoms,
        size_t neval)
{
    std::cout << std::endl; // flush the output
    std::cout << "======================================================\n";
    std::cout << "benchmarking " << neval << " energy evaluations for " << natoms << " atoms\n";


    pele::Array<double> x(3*natoms);

    Timer timer;
    timer.start();
    auto pot = pot_maker.get_potential_coords(x);
    timer.stop();
    double time_init = timer.get();
    std::cout <<  "  time to initialize potential: " << time_init << "\n";

    double t = bench_potential(pot, x, neval);
    std::cout << "  total time: " << t << "\n";
    std::cout << "  timer per evaluation: " << t / neval << "\n";
    std::cout << "  timer per evaluation per atom: " << t / neval / natoms << "\n";
    std::cout << "  natoms neval time_init time_eval\n";
    std::cout << "  all times:"
            << " " << natoms
            << " " << neval
            << " " << time_init
            << " " << t
            << "\n";
    return t;
}
