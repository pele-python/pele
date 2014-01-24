#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <ctime>
#include "base_potential.h"
#include "array.h"
#include "optimizer.h"
#include "base_integrator.h"
#include "velocity_verlet.h"
#include "lj.h"
#include "_modified_fire.h"

int main(int argc, char ** argv)
{
	double r;
	pele::Array<double> x(9);
	pele::Array<double> grad(9);
	pele::Array<double> v(9);
	srand(time(NULL));

	for(int i=0; i<x.size();++i)
	{
		r = ((double) rand() / (RAND_MAX)) + 1;
		x[i] = r;
	}

	pele::LJ potential(1,1);
	pele::BasePotential * potentialptr;

	potentialptr = &potential;

	cout<<"potentialptr address: "<<potentialptr<<std::endl;

	potentialptr->get_energy_gradient(x, grad);

	for(int i=0; i<x.size();++i)
		{
			cout<<"grad "<<grad[i]<<"\n";
		}

	double dt = 0.1;

	assert(!x.empty());

	pele::VelocityVerlet integrator(potentialptr, x, dt);
	integrator.wrap_v(v);

	std::cout<<v<<std::endl;
	integrator.oneiteration();
	std::cout<<v<<std::endl;
	std::cout<<"end of integrator integrator test"<<std::endl;

	pele::MODIFIED_FIRE minimiser(potentialptr, x, 0.1, 1);

	minimiser.run(100000000);
	//std::cout<<v<<std::endl;


	return 0;
}
