#include <iostream>
#include <assert.h>
#include "neb.h"

namespace pele {

void NEB::initialize(void)
{
	_nimages = 30;
}

void NEB::set_nimages(int nimages)
{
	_nimages = nimages;
}

void NEB::set_path(Array &x1, Array &x2)
{
	std::list< Array * > path;
	path.push_back(&x1);
	path.push_back(&x2);

	set_path(path);
}

void NEB::set_path(std::list< Array * > &path)
{
	double s=0;
	auto iter1 = path.begin();
	auto iter2 = path.begin();
	int n = (*iter1)->size();

	for(iter2++; iter2!=path.end(); ++iter2) {
		s += _neb_distance->dist(**iter1, **iter2);
		iter1 = iter2;
	}
	std::cout << "The total path distance is:" << s << std::endl;

	_neb_coords.resize(n*_nimages);
}

void NEB::interpolate(Array &x1, Array &x2, Array &xout, double t)
{
	assert(x1.size() == x2.size());
	assert(x1.size() == xout.size());

	for(int i=0;i<x1.size();++i) {
		xout[i] = t*x2[i] - (1.-t)*x1[i];
	}
}
void test_array(Array a)
{
	std::cout << a << std::endl;
}

void test_py(void)
{
	std::cout << "hello world\n";
}
}
