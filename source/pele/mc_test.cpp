#include <iostream>
#include "harmonic.h"
#include "mc.h"
#include "actions.h"
#include "histogram.h"
#include "takestep.h"
#include "accept_test.h"
#include <memory>

/*returns the energy of an harmonic oscillator*/

using std::shared_ptr;

int main(){

  int ntot =1000000;
  pele::Array<double> coords(3,0);
  pele::Array<double> origin(3,0);
  pele::Histogram * hist = new pele::Histogram(2,15,0.1);
  pele::RandomCoordsDisplacement * random = new pele::RandomCoordsDisplacement(coords.size());
  pele::MetropolisTest * metropolis = new pele::MetropolisTest;
  pele::RecordEnergyHistogram * histogram = new pele::RecordEnergyHistogram(hist);
  pele::Harmonic * harmonic = new pele::Harmonic(origin, 1);
  pele::MC MC(harmonic, coords, 1, 1);
  MC.add_action(shared_ptr<pele::RecordEnergyHistogram>(histogram));
  MC.add_accept_test(shared_ptr<pele::MetropolisTest>(metropolis));
  MC.set_takestep(shared_ptr<pele::RandomCoordsDisplacement>(random));
  std::cout<<"this is a test"<<std::endl;
  MC.run(ntot);
  hist->print(ntot);
  std::vector<size_t>::iterator it;
  size_t renorm = 0;
  for(it = hist->begin();it != hist->end();++it)
  {
	  renorm += *it;
  }
  std::cout<<"renorm "<<renorm<<std::endl;
  std::cout<<"histogram niter "<<hist->_niter<<std::endl;
  renorm = 0;
  /*for(it = (*hist).begin();it != (*hist).end();++it)
    {
  	  *it /= (renorm*exp(-*it));
  	  //std::cout<<*it<<std::endl;
    }*/
  /*delete random;
  delete metropolis;
  delete histogram;
  delete harmonic;*/

  return 0;

}
