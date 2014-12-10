#include <random>
#include <iostream>
#include <fstream>
#include <string>

#include "pele/lj.h"
#include "pele/lbfgs.h"
#include "pele/matrix.h"
#include "pele/array.h"

#include "bench_utils.hpp"

using namespace pele;
using std::string;

int main(int argc, char ** argv)
{
    std::cout << std::setprecision(16);
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0,1);

    size_t natoms = 200;
    Array<double> x(3*natoms);



    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = 2 * (1 - 0.5 * distribution(generator));
    }

    if (argc > 1) {
        auto fname = string(argv[1]);
        std::cout << "reading from file " << fname << std::endl;
        x = coords_from_file(fname);
        natoms = x.size() / 3;
    }

    // print the coords array array
//    pele::MatrixAdapter<double> m(x, 3);
//    std::cout << "initial coordinates\n";
//    std::cout << m << std::endl;

    auto lj = std::make_shared<LJ>(4., 4.);

    std::cout << "energy " << lj->get_energy(x) << std::endl;

    LBFGS lbfgs(lj, x);
    lbfgs.set_max_iter(100000);
    lbfgs.set_iprint(1000);
    lbfgs.run();

}
