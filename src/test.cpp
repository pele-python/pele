#include <vector>
#include <iostream>
#include <math.h>
#include "neb.h"

using namespace pygmin;

class Dist : public NEB::NEBDistance
{
public:
    double dist(Array &x1, Array &x2) {
    	double d=0;
    	for(int i=0; i<x1.size(); ++i) {
    		double tmp = x2[i] - x1[i];
    		d+=tmp*tmp;
    	}
    	return sqrt(d);
    }
};

class MyPotential : public Potential
{
public:
	double get_energy(Array &x) { return 0.; }
	double get_energy_gradient(Array &x, Array &grad) {
		grad = 0.;
		return 0.;
	}
};

int main(int argc, char **argv)
{
	Dist dist;
	MyPotential pot;

	std::vector<double> x1 = { 1., 0. };
	std::vector<double> x2 = { 0., 1. };

	NEB neb(&pot, &dist);
	neb.set_nimages(10);

	Array ax1(x1);
	Array ax2(x2);

	neb.set_path(ax1, ax2);
}
