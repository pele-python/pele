#ifndef _pele_vecn_h_
#define _pele_vecn_h_

#include <cmath>
#include <stdlib.h>

namespace pele{
template<size_t N>
class VecN {
    typedef double dtype;
    dtype m_data[N];

public:

    /**
     * default constructor
     */
    VecN(){}

    /**
     * initialize with constant
     */
    VecN(dtype const & d) { assign(d); }

    size_t size() const { return N; }

    /**
     * return pointer to data
     */
    inline dtype * data() { return m_data; }
    inline dtype const * data() const { return m_data; }

    /**
     * return iterators over data
     */
    typedef dtype * iterator;
    typedef dtype const * const_iterator;
    inline iterator begin() { return data(); }
    inline iterator end() { return data() + size(); }
    inline const_iterator begin() const { return data(); }
    inline const_iterator end() const { return data() + size(); }

    /**
     * access an element in the vector
     */
    inline dtype & operator[](const size_t i) { return m_data[i]; }
    inline dtype const & operator[](const size_t i) const { return m_data[i]; }

    /**
     * assign each element of the vector to be d
     */
    void assign(dtype const & d)
    {
        for (size_t i=0; i<N; ++i){
            m_data[i] = d;
        }
    }


    /*
     * Compound Assignment Operators += -= *=
     */
    VecN<N> &operator+=(const VecN<N> & rhs) {
        for (size_t i = 0; i < N; ++i) {
            m_data[i] += rhs[i];
        }
        return *this;
    }

    VecN<N> &operator+=(const dtype &rhs) {
        for (size_t i = 0; i < N; ++i) {
            m_data[i] += rhs;
        }
        return *this;
    }

    VecN<N> &operator-=(const VecN<N> & rhs){
        for (size_t i = 0; i < N; ++i) {
            m_data[i] -= rhs[i];
        }
        return *this;
    }

    VecN<N> &operator-=(const dtype &rhs) {
        for (size_t i = 0; i < N; ++i) {
            m_data[i] -= rhs;
        }
        return *this;
   }

    VecN<N> &operator*=(const VecN<N> & rhs){
        for (size_t i = 0; i < N; ++i) {
            m_data[i] *= rhs[i];
        }
        return *this;
    }

    VecN<N> &operator*=(const dtype &rhs) {
        for (size_t i = 0; i < N; ++i) {
            m_data[i] *= rhs;
        }
        return *this;
    }


    VecN<N> &operator/=(const VecN<N> & rhs){
        for (size_t i = 0; i < N; ++i) {
            m_data[i] /= rhs[i];
        }
        return *this;
    }

    VecN<N> &operator/=(const  dtype &rhs) {
        for (size_t i = 0; i < N; ++i) {
            m_data[i] /= rhs;
        }
        return *this;
    }


    /**
     * returns the sum of all elements (reduces the array)
     */
    const dtype sum() const {
        dtype sum_array = 0;
        for (size_t i = 0; i<N; ++i){
            sum_array += m_data[i];
        }
        return sum_array;
    }

    /**
     * returns the product of all elements (reduces the array)
     */
    const dtype prod() const {
        dtype p = 1;
        for (size_t i = 0; i<N; ++i){
            p *= m_data[i];
        }
        return p;
    }

};

/**
 * compute the dot product of two Arrays
 */
template<size_t N>
inline double dot(VecN<N> const &v1, VecN<N> const &v2)
{
  double dot = 0.;
  for (size_t i=0; i<N; ++i) {
    dot += v1[i] * v2[i];
  }
  return dot;
}

/**
 * compute the L2 norm of an Array
 */
template<size_t N>
inline double norm(VecN<N> const &v)
{
  return sqrt(dot(v, v));
}


} // close namespace pele
#endif
