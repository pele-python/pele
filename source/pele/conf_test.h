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
	pele::periodic_distance _periodic_dist;
	pele::cartesian_distance _cartesian_dist;
	pele::GradientOptimizer * _optimizer;
	Array<double> _origin, _hs_radii, _rattlers, _distance;
	double _dtol, _d, _rms;
	size_t _N, _Nnoratt, _nparticles, _ialign;
public:
	CheckSameMinimum(pele::GradientOptimizer * optimizer, Array<double> origin, Array<double> hs_radii, Array<double> boxvec,
			Array<double> rattlers, double dtol);
	virtual bool test(Array<double> &trial_coords, MC * mc);
	virtual ~CheckSameMinimum(){
		if (_optimizer != NULL)
			delete _optimizer;
	}
	double get_distance(){return _d;}
	Array<double> get_distance_array(){
		Array<double> x(_distance.copy());
		return x;
	}
	inline bool check_overlap(Array<double> &trial_coords);
	inline void align(Array<double> &trial_coords);
	inline void get_periodic_distance(Array<double> &trial_coords);
};

CheckSameMinimum::CheckSameMinimum(pele::GradientOptimizer * optimizer, Array<double> origin, Array<double> hs_radii,
		Array<double> boxvec, Array<double> rattlers, double dtol):
		_periodic_dist(boxvec[0], boxvec[1], boxvec[2]),
		_cartesian_dist(pele::cartesian_distance()),
		_optimizer(optimizer), _origin(origin), _hs_radii(hs_radii.copy()),
		_rattlers(rattlers), _distance(origin.size(),0),_dtol(dtol),_d(0),
		_rms(0),_N(origin.size()),_Nnoratt(0), _nparticles(_N/3), _ialign(0){
			for(size_t i=0;i<_N;++i){_Nnoratt += _rattlers[i];}
			//also find first non rattler
			for(size_t i=0;i<_N;++i){
				if (_rattlers[i] == 1){
					_ialign=3; //i ->3 is only for testing reasons!!!!!!!!!
					break;}
			}
		}

inline void CheckSameMinimum::align(Array<double> &trial_coords){
	size_t i,j, i1;
	double dr[3];

	//it shouldn't matter whether one uses periodic or cartesian distance here
	//because the rms displacement at the end of the test is che
	_cartesian_dist.get_rij(dr, &_origin[_ialign*3], &trial_coords[_ialign*3]);
	//std::cout<<"dr "<<dr[0]<<" "<<dr[1]<<" "<<dr[2]<<std::endl;
	for(i=0;i<_nparticles;++i){
			i1 = 3*i;
			for(j=0;j<3;++j){
				trial_coords[i1+j] += dr[j];
			}
	}
}

inline void CheckSameMinimum::get_periodic_distance(Array<double> &trial_coords){
	size_t i,i1;
	double dr[3];

	for(i=0;i<_nparticles;++i){
		i1 = 3*i;
		_periodic_dist.get_rij(dr, &_origin[i1], &trial_coords[i1]);
		_distance[i1] = dr[0];
		_distance[i1+1] = dr[1];
		_distance[i1+2] = dr[2];
		}
}


inline bool CheckSameMinimum::check_overlap(Array<double> &trial_coords){
	size_t i,j, i1, j1;
	double dr[3];
	double dij2,dij;

	for(i=0;i<_nparticles;++i){
		i1 = 3*i;
		for(j=0;j<_nparticles;++j){
			if(i != j){
				j1 = 3*j;
				_periodic_dist.get_rij(dr, &trial_coords[i1], &trial_coords[j1]);
				dij2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
				dij = sqrt(dij2);
				dij -= (_hs_radii[i] + _hs_radii[j]);
				if (dij <= 0)
					return false;
				}
			}
		}

	return true;
}

bool CheckSameMinimum::test(Array<double> &trial_coords, MC * mc)
{
	bool quench_success, no_overlap;
	size_t nfev;

	no_overlap = this->check_overlap(trial_coords);

	if (! no_overlap)
		return false;

	_optimizer->reset(trial_coords);
	_optimizer->run();

	//add number of energy evaluations to mc eval count
	nfev = _optimizer->get_niter();
	mc->_neval += nfev;

	//first test: minimisation must have converged
	quench_success = _optimizer->success();
	if (! quench_success)
		return false;

	//copy coordinates of quenched structure into _distance array
	_distance.assign(_optimizer->get_x());
	//align
	this->align(_distance);

	//compute distances subtracting the origin's coordinates using PBC
	this->get_periodic_distance(_distance);

	//set to 0 distances of rattlers
	for (size_t l = 0; l < _N; ++l){
		_distance[l] *= _rattlers[l];
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
