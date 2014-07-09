#ifndef _PELE_DISTANCE_H
#define _PELE_DISTANCE_H

#include <cmath>
#include <stdexcept>
#include "array.h"

// roundis missing in visual studio
#ifdef _MSC_VER
    inline double round(double r) {
        return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
    }
#endif

/**
 * These classes and structs are used by the potentials to compute distances.
 * They must have a member function get_rij() with signature
 *
 *    void get_rij(double * r_ij, double const * const r1, 
 *                                double const * const r2) 
 *
 * Where r1 and r2 are the position of the two atoms and r_ij is an array of
 * size 3 which will be used to return the distance vector from r1 to r2.
 */

namespace pele
{

/**
* compute the cartesian distance
*/

template<size_t ndim>
struct cartesian_distance {
    static const size_t _ndim = ndim;
    inline void get_rij(double * const r_ij, double const * const r1,
            double const * const r2) const
    {
        for (size_t i=0;i<ndim;++i)
            r_ij[i] = r1[i] - r2[i];
    }

};

//cartesian distance template specializations

template<>
struct cartesian_distance <3> {
    static const size_t _ndim = 3;
    inline void get_rij(double * const r_ij, double const * const r1,
                 double const * const r2) const
    {
        r_ij[0] = r1[0] - r2[0];
        r_ij[1] = r1[1] - r2[1];
        r_ij[2] = r1[2] - r2[2];
    }
};

template<>
struct cartesian_distance <2> {
    static const size_t _ndim = 2;
    inline void get_rij(double * const r_ij, double const * const r1,
               double const * const r2) const
    {
        r_ij[0] = r1[0] - r2[0];
        r_ij[1] = r1[1] - r2[1];
    }
};

/**
* periodic boundary conditions in rectangular box
*/

template<size_t ndim>
class periodic_distance {
public:
    static const size_t _ndim = ndim;
    double _box[ndim];
    double _ibox[ndim];

    periodic_distance(Array<double> const box)
    {
        if (box.size() != _ndim) {
            throw std::invalid_argument("box.size() must be equal to ndim");
        }
        for (size_t i=0;i<ndim;++i) {
            _box[i] = box[i];
            _ibox[i] = 1/box[i];
        }
    }

    periodic_distance()
    { 
        throw std::runtime_error("the empty constructor is not available for periodic boundaries"); 
    }

    inline void get_rij(double * const r_ij, double const * const r1,
                 double const * const r2) const
    {
        for (size_t i=0;i<ndim;++i) {
            r_ij[i] = r1[i] - r2[i];
            r_ij[i] -= round(r_ij[i] * _ibox[i]) * _box[i];
        }
    }
};

//periodic distance template specializations

template <>
class periodic_distance <3> {
public:
    static const size_t _ndim = 3;
    double const _boxx;
    double const _boxy;
    double const _boxz;
    double const _iboxx;
    double const _iboxy;
    double const _iboxz;

    periodic_distance(Array<double> const box)
        : _boxx(box[0]),
          _boxy(box[1]),
          _boxz(box[2]),
          _iboxx(1./_boxx),
          _iboxy(1./_boxy),
          _iboxz(1./_boxz)
    {
        if (box.size() != _ndim) {
            throw std::invalid_argument("box.size() must be equal to ndim");
        }
    }

    /* this constructor exists only so the compiler doesn't complain.
    * It should never be used */
    periodic_distance()
        : _boxx(), _boxy(), _boxz(), _iboxx(), _iboxy(), _iboxz()
    { 
        throw std::runtime_error("the empty constructor is not available for periodic boundaries"); 
    }

    inline void get_rij(double * const r_ij, double const * const r1, 
                 double const * const r2) const
    {
        r_ij[0] = r1[0] - r2[0];
        r_ij[1] = r1[1] - r2[1];
        r_ij[2] = r1[2] - r2[2];
        r_ij[0] -= round(r_ij[0] * _iboxx) * _boxx;
        r_ij[1] -= round(r_ij[1] * _iboxy) * _boxy;
        r_ij[2] -= round(r_ij[2] * _iboxz) * _boxz;
    }
};

template <>
class periodic_distance <2>{
public:
    static const size_t _ndim = 2;
    double const _boxx;
    double const _boxy;
    double const _iboxx;
    double const _iboxy;

    periodic_distance(Array<double> const box)
        : _boxx(box[0]),
          _boxy(box[1]),
          _iboxx(1./_boxx),
          _iboxy(1./_boxy)
    {
        if (box.size() != _ndim) {
            throw std::invalid_argument("box.size() must be equal to ndim");
        }
    }

    /* this constructor exists only so the compiler doesn't complain.
     * It should never be used */
    periodic_distance()
        : _boxx(), _boxy(), _iboxx(), _iboxy()
    { throw std::runtime_error("the empty constructor is not available for periodic boundaries"); }

    inline void get_rij(double * const r_ij, double const * const r1,
                double const * const r2) const
    {
        r_ij[0] = r1[0] - r2[0];
        r_ij[1] = r1[1] - r2[1];
        r_ij[0] -= round(r_ij[0] * _iboxx) * _boxx;
        r_ij[1] -= round(r_ij[1] * _iboxy) * _boxy;
    }
};

/*
 * Interface classes to pass a generic non template pointer to DistanceInterface
 * (this should reduce the need for templating where best performance is not essential)
 */

class DistanceInterface{
protected:
public:
    virtual void get_rij(double * const r_ij, double const * const r1, double const * const r2) const =0;
    virtual ~DistanceInterface(){ }
};

template<size_t ndim>
class CartesianDistanceWrapper : public DistanceInterface{
protected:
    cartesian_distance<ndim> _dist;
public:
    static const size_t _ndim = ndim;
    inline void get_rij(double * const r_ij, double const * const r1,
            double const * const r2) const
    {
        _dist.get_rij(r_ij, r1, r2);
    }
};

template<size_t ndim>
class PeriodicDistanceWrapper : public DistanceInterface{
protected:
    periodic_distance<ndim> _dist;
public:
    static const size_t _ndim = ndim;
    PeriodicDistanceWrapper(Array<double> const box)
        : _dist(box)
    {};

    inline void get_rij(double * const r_ij, double const * const r1,
            double const * const r2) const
    {
        _dist.get_rij(r_ij, r1, r2);
    }
};

}
#endif
