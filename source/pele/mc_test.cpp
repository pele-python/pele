#include <iostream>
#include "harmonic.h"
#include "mc.h"
#include "actions.h"
#include "histogram.h"
#include "takestep.h"
#include "accept_test.h"
#include <fstream>
#include <memory>

/*returns the energy of an harmonic oscillator*/

using std::shared_ptr;

int main(){

  double Emin = 0.5;
  double Emax = 2;
  int ntot =10000000;
  double bin = 0.005;
  double temperature = 0.8;
  double stepsize = 1;
  int ndim = 5;
  pele::Array<double> coords(ndim,0);
  pele::Array<double> origin(ndim,0);
  pele::Histogram * hist = new pele::Histogram(0,25,bin);

  pele::RandomCoordsDisplacement * random = new pele::RandomCoordsDisplacement(coords.size());
  pele::MetropolisTest * metropolis = new pele::MetropolisTest;
  pele::RecordEnergyHistogram * histogram = new pele::RecordEnergyHistogram(hist);
  pele::Harmonic * harmonic = new pele::Harmonic(origin, 1);
  pele::MC MC(harmonic, coords, temperature, stepsize);
  MC.add_action(shared_ptr<pele::RecordEnergyHistogram>(histogram));
  MC.add_accept_test(shared_ptr<pele::MetropolisTest>(metropolis));
  MC.set_takestep(shared_ptr<pele::RandomCoordsDisplacement>(random));

  //generate initial coordinates
  pele::RandomCoordsDisplacement random_init = pele::RandomCoordsDisplacement(coords.size());
  random_init.takestep(coords, sqrt(2*Emax)*ndim/(ndim+1), &MC); //expectation value for a point in harmonic well

  std::cout<<"this is a test"<<std::endl;
  MC.run(ntot);
  std::cout<<"HISTOGRAM"<<std::endl;
  //hist->print(ntot);
  std::vector<size_t>::iterator it;
  std::cout<<"histogram niter "<<hist->_niter<<std::endl;
  std::cout<<"accepted fraction "<<MC.get_accepted_fraction()<<std::endl;

  std::ofstream myfile;
  myfile.open("harmonic_dos.dat");

  double h;
  double E = 0;
  double sumh = 0;

  for(it = hist->begin();it != hist->end();++it)
  {
	  if (E >= Emin && E <= Emax)
	  {
	  sumh += ((double) *it) * exp(temperature*E) * bin;
	  }
	  E += bin;
  }

  E = 0;

  for(it = hist->begin();it != hist->end();++it)
    {
  	  h = ((double) *it) * exp(temperature*E) / sumh;
  	  myfile<<E<<"\t"<<h<<std::endl;
  	  E += bin;
  	  //myfile<<E<<"\t"<<h<<"\t"<<std::endl;
    }

  myfile.close();

  myfile.open("analytical_dos.dat");

  E = 0;
  sumh = 0;

    for(it = hist->begin();it != hist->end();++it)
    {
    	if (E >= Emin && E <= Emax)
    	{
    	sumh += pow(E,ndim/2-1) * bin;
    	}
    	E += bin;
    }

  E = 0;

  for(it = hist->begin();it != hist->end();++it)
  {
	  h = pow(E,ndim/2-1)/sumh;
	  myfile<<E<<"\t"<<h<<std::endl;
	  E += bin;
  }

  myfile.close();
  return 0;

}
