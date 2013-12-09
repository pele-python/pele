#ifndef _PELE_DISTANCE_H
#define _PELE_DISTANCE_H

#include <cmath>
#include <stdexcept>

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
  struct cartesian_distance {
    inline void get_rij(double *  r_ij, 
                 double const * const r1, 
                 double const * const r2) 
    {
      r_ij[0] = r1[0] - r2[0];
      r_ij[1] = r1[1] - r2[1];
      r_ij[2] = r1[2] - r2[2];
    } 
  };

  /**
   * periodic boundary conditions in rectangular box
   */
  class periodic_distance {
    public: 
      double const _boxx;
      double const _boxy;
      double const _boxz;
      double const _iboxx;
      double const _iboxy;
      double const _iboxz;

      periodic_distance(double boxx, double boxy, double boxz) :
          _boxx(boxx),
          _boxy(boxy),
          _boxz(boxz),
          _iboxx(1./_boxx),
          _iboxy(1./_boxy),
          _iboxz(1./_boxz)
      {}

      /* this constructor exists only so the compiler doesn't complain.
       * It should never be used */
      periodic_distance() :
          _boxx(1000.),
          _boxy(1000.),
          _boxz(1000.),
          _iboxx(1./_boxx),
          _iboxy(1./_boxy),
          _iboxz(1./_boxz)
      { throw std::runtime_error("the empty constructor is not available for periodic boundaries"); }

      inline void get_rij(double * r_ij, 
                     double const * const r1, 
                     double const * const r2) 
      {
          r_ij[0] = r1[0] - r2[0];
          r_ij[1] = r1[1] - r2[1];
          r_ij[2] = r1[2] - r2[2];
          r_ij[0] -= round(r_ij[0] * _iboxx) * _boxx;
          r_ij[1] -= round(r_ij[1] * _iboxy) * _boxy;
          r_ij[2] -= round(r_ij[2] * _iboxz) * _boxz;
      } 
  };
}
#endif
