#ifndef _PELE_CONF_TEST_H__
#define _PELE_CONF_TEST_H__

#include <iostream>
#include <math.h>
#include <algorithm>
#include <random>
#include <chrono>
#include "array.h"
#include "mc.h"
#include "optimizer.h"
#include "distance.h"

using std::runtime_error;
using pele::Array;

namespace pele{

class CheckSphericalContainer:public ConfTest{
protected:
	double _radius2;
public:
	CheckSphericalContainer(double radius);
	virtual bool test(Array<double> &trial_coords, MC * mc);
	virtual ~CheckSphericalContainer(){}
};

CheckSphericalContainer::CheckSphericalContainer(double radius):
		_radius2(radius*radius){}

bool CheckSphericalContainer::test(Array<double> &trial_coords, MC * mc)
{
  double r2;
  int i, j, N;
  N = trial_coords.size();

  for (i=0; i<N; i+=3)
  {

	  r2 = 0;
	  for (j=i; j<i+3; ++j)
	  {
		r2 += trial_coords[j] * trial_coords[j];
	  }
	  if (r2 > _radius2)
	  {
	    //printf("fail spherical container %d %f %f %f %f\n", i, sqrt(r2), x[i], x[i+1], x[i+2]);
	    //an atom is outside the spherical container
	    return false;
	  }
  }

  //printf("check spherical OK ");
  return true;
}

}


#endif
