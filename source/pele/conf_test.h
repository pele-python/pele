#ifndef _PELE_CONF_TEST_H__
#define _PELE_CONF_TEST_H__

#include <math.h>
#include <algorithm>
#include <random>
#include <chrono>
#include "array.h"
#include "mc.h"
#include "optimizer.h"

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


/*check same minimum class
 * _optimizer: pointer to object of class GradientOptimizer performing minimisation according to some potential
 * 				passed to the object during its construction
 * _origin: coordinates to which the quenched structure is compared to
 * _rattlers: array of 1s or 0s: if not indicates a rattler,0 -> rattler
 * 															1 -> jammed particle
 * 			this convention removes if statements in the for loop and replaces them with
 * 			arithmetic operation (distance[i] *= rattlers[i]), distance is set artificially to
 * 			zero if the particle is a rattler. note _rattlers.size() = coords.size()
 * _distance: array containing the Euclidean distance between trial_coords and origin
 * _d: norm of distance
 * _rms: root mean square displacement from origin
 * _E = energy of the quenched state
 * _dtol: tolerance on distances
 * _Etol: tolerance on energies (a minimum should be whithin this value from _Emin)
 * _Eor: energy of the origin (must pass it because CheckSameMinimum knows nothing about the potential used by the optimiser)
 * _Nnoratt: total number of non rattlers degrees of freedom
 * */

class CheckSameMinimum:public ConfTest{
protected:
	pele::GradientOptimizer * _optimizer;
	Array<double> _origin, _rattlers, _distance;
	double _Eor, _Etol, _dtol, _E, _d, _rms;
	int _Nnoratt;
public:
	CheckSameMinimum(pele::GradientOptimizer * optimizer, Array<double> origin, Array<double> rattlers, double Eor, double Etol, double dtol);
	virtual bool test(Array<double> &trial_coords, MC * mc);
	virtual ~CheckSameMinimum(){}
	double get_distance(){return _d;}
	Array<double> get_distance_array(){
		Array<double> x(_distance.copy());
		return x;
	}
};

CheckSameMinimum::CheckSameMinimum(pele::GradientOptimizer * optimizer, Array<double> origin, Array<double> rattlers,
		double Eor, double Etol, double dtol):
		_optimizer(optimizer), _origin(origin), _rattlers(rattlers),
		_distance(origin.size()),_Eor(Eor),_Etol(Etol),_dtol(dtol),
		_E(_Eor),_d(0),_rms(0),_Nnoratt(0){
		for(int i=0;i<_rattlers.size();++i){
			_Nnoratt += _rattlers[i];
		}
}

bool CheckSameMinimum::test(Array<double> &trial_coords, MC * mc)
{
	bool quench_success;

	_optimizer->reset(trial_coords);
	_optimizer->run();
	_E = _optimizer->get_f();
	quench_success = _optimizer->success();
	//first test: minimisation must have converged and energy must be within some reasonable range of Eor
	if ((quench_success == false) || (abs(_E - _Eor) > _Etol))
	  return false;

	//copy coordinates of quenched structure
	_distance = (_optimizer->get_x()).copy();
	//compute distances subtracting the origin's coordinates
	_distance -= _origin;
	//set to 0 distances of rattlers
	size_t N = _distance.size();
	for (size_t j = 0; j < N; ++j){
		_distance[j] *= _rattlers[j];
	}
	//compute rms displacement from origin
	_d = norm(_distance);
	_rms = _d / sqrt(_Nnoratt);
	if (_rms > _dtol)
		return false;
	else
		return true;
}

}


#endif
