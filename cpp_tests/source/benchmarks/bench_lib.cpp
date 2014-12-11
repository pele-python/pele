#include <random>
#include <iostream>
#include <fstream>
#include <string>

#include "pele/neighbor_iterator.h"
#include "pele/lj_cut.h"
#include "pele/lbfgs.h"
#include "pele/matrix.h"


typedef pele::LJCutPeriodicCellLists<3> LJCutPeriodicCellLists3;
typedef pele::Array<double> ArrayD;
typedef pele::LBFGS LBFGS_def;
