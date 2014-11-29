#include <random>
#include <iostream>

#include "pele/neighbor_iterator.h"
#include "pele/lj_cut.h"
#include "pele/lbfgs.h"

using namespace pele;

int main()
{
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0,1);

    size_t natoms = 1600;
    double rcut = 2.;
    double density = 1.2;
    double boxl = std::pow(natoms / density * (4/3 * M_PI), 1./3);
    std::cout << "box length " << boxl << std::endl;
    Array<double> boxvec(3, boxl);

    Array<double> x(3*natoms);
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = distribution(generator) * boxl;
    }

    double ncellx_scale = 1.;
    auto lj = std::make_shared<LJCutPeriodicCellLists<3> >(4., 4., rcut, boxvec, ncellx_scale);

    std::cout << "energy " << lj->get_energy(x) << "\n";

    LBFGS lbfgs(lj, x);
    lbfgs.set_max_iter(100000);
    lbfgs.set_iprint(100);
    lbfgs.run();

}
