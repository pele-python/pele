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

  pele::Array<double> coords(3,0);
  pele::Array<double> origin(3,0);
  pele::Histogram * hist = new pele::Histogram(2,9,0.1);
  pele::RandomCoordsDisplacement * random = new pele::RandomCoordsDisplacement(coords.size());
  pele::MetropolisTest * metropolis = new pele::MetropolisTest;
  pele::RecordEnergyHistogram * histogram = new pele::RecordEnergyHistogram(hist);
  pele::Harmonic * harmonic = new pele::Harmonic(origin, 1);
  pele::MC MC(harmonic, coords, 1, 1);
  MC.add_action(shared_ptr<pele::RecordEnergyHistogram>(histogram));
  MC.add_accept_test(shared_ptr<pele::MetropolisTest>(metropolis));
  MC.set_takestep(shared_ptr<pele::RandomCoordsDisplacement>(random));
  std::cout<<"this is a test"<<std::endl;
  MC.run(1000);
  hist->print();


return 0;
}
