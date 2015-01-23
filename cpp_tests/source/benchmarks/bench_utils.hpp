#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <ctime>


#include "pele/array.h"

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
