#ifndef _pele_vecn_h_
#define _pele_vecn_h_

#include <cmath>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <pele/array.h>
#include <string>
#include <sstream>
#include <algorithm>

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

    /**
     * initialize as copy of pele array
     */
    VecN(pele::Array<dtype> const & x)
    {
        if (x.size() != N) {
            throw std::runtime_error("VecN constructor: array must have the same size as vector");
        }
        for (size_t i = 0; i < N; ++i) {
            m_data[i] = x[i];
        }
    }

    /**
     * initialize values from input iterators
     */
    template<class input_iter>
    VecN(input_iter ibegin, input_iter iend) {
        std::copy(ibegin, iend, this->begin());
    }

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

    /**
     * copy the data in a pele::Array into this vector
     */
    VecN<N> & operator=(pele::Array<dtype> const & rhs) {
        if (rhs.size() != N) {
            throw std::runtime_error("operator=: array must have the same size");
        }
        for (size_t i = 0; i < N; ++i) {
            m_data[i] = rhs[i];
        }
        return *this;
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

    VecN<N> operator-(VecN<N> const & rhs) const {
        VecN<3> v;
        for (size_t i = 0; i < N; ++i) {
            v[i] = m_data[i] - rhs[i];
        }
        return v;
    }


    /**
     * returns the sum of all elements (reduces the array)
     */
    dtype sum() const {
        dtype sum_array = 0;
        for (size_t i = 0; i<N; ++i){
            sum_array += m_data[i];
        }
        return sum_array;
    }

    /**
     * returns the product of all elements (reduces the array)
     */
    dtype prod() const {
        dtype p = 1;
        for (size_t i = 0; i<N; ++i){
            p *= m_data[i];
        }
        return p;
    }

}; // close VecN

template<size_t N, size_t M>
class MatrixNM {
    typedef double dtype;
    static size_t const m_size = N * M;
    dtype m_data[m_size];

public:

    /**
     * default constructor
     */
    MatrixNM() {}

    /**
     * initialize with constant
     */
    MatrixNM(dtype const & d) { assign(d); }

    /**
     * initialize as copy of pele array
     */
    MatrixNM(pele::Array<dtype> const & x)
    {
        if (x.size() != m_size) {
            std::stringstream ss;
            ss << "MatrixNM constructor: array (size " << x.size() << ") must have the same size as matrix " << m_size;
            throw std::runtime_error(
                    ss.str()
                    );
        }
        for (size_t i = 0; i < m_size; ++i) {
            m_data[i] = x[i];
        }
    }

    size_t size() const { return m_size; }

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
     * assign each element of the vector to be d
     */
    void assign(dtype const & d)
    {
        for (size_t i=0; i<m_size; ++i){
            m_data[i] = d;
        }
    }

    /**
     * provide access to matrix element at row i and column j
     */
    inline dtype const & operator()(size_t i, size_t j) const
    {
        return m_data[i * M + j];
    }
    inline dtype & operator()(size_t i, size_t j)
    {
        return m_data[i * M + j];
    }

    inline std::pair<size_t, size_t> shape() const
    {
        return std::pair<size_t, size_t>(N, M);
    }

    MatrixNM<N, M> &operator*=(dtype const & rhs) {
        for (size_t i = 0; i < m_size; ++i) {
            m_data[i] *= rhs;
        }
        return *this;
    }

    dtype trace()
    {
        dtype t = 0;
        for (size_t i = 0; i<N; ++i){
            t += (*this)(i,i);
        }
        return t;
    }

    MatrixNM<N,M> operator-(MatrixNM<N,M> const & rhs) const {
        MatrixNM<N,M> v;
        for (size_t i = 0; i < m_size; ++i) {
            v.m_data[i] = m_data[i] - rhs.m_data[i];
        }
        return v;
    }


}; // close MatrixNM

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


/**
 * matrix_multiplication
 *
 * This is a really simple implementation of matrix multiplication.
 * It is order N*M*L but can be done with much better scaling
 */
template<size_t N, size_t L, size_t M>
MatrixNM<N,M> dot(MatrixNM<N,L> const & A, MatrixNM<L, M> const & B)
{
    MatrixNM<N,M> C(0);
    for (size_t i = 0; i<N; ++i){
        for (size_t j = 0; j<M; ++j){
            double val = 0;
            for (size_t k = 0; k<L; ++k){
                val += A(i,k) * B(k,j);
            }
            C(i,j) = val;
        }
    }
    return C;
}

/**
 * multiply a matrix times an vector
 */
template<size_t N, size_t M>
pele::VecN<N> dot(MatrixNM<N,M> const & A, pele::VecN<M> const & v)
{
    pele::VecN<N> C(0);
    for (size_t i = 0; i<N; ++i){
        double val = 0;
        for (size_t k = 0; k<M; ++k){
            val += A(i,k) * v[k];
        }
        C[i] = val;
    }
    return C;
}

template<size_t N, size_t M>
pele::MatrixNM<M,N> transpose(MatrixNM<N,M> const & A)
{
    pele::MatrixNM<M,N> mat;
    for (size_t i = 0; i<N; ++i){
        for (size_t k = 0; k<M; ++k){
            mat(k,i) = A(i,k);
        }
    }
    return mat;
}

template<size_t N>
pele::MatrixNM<N,N> identity()
{
    pele::MatrixNM<N,N> A(0.);
    for (size_t i = 0; i<N; ++i) {
        A(i,i) = 1.;
    }
    return A;
}

// for matrix printing
template<size_t N, size_t M>
inline std::ostream &operator<<(std::ostream &out, const MatrixNM<N,M> &a) {
    out << "[ ";
    for(size_t n=0; n<N;++n) {
        for(size_t m=0; m<M;++m) {
            if(m>0) out << ", ";
            out << a(n,m);
        }
        if (n < N-1) out << ",\n  ";
    }
    out << " ]";
    return out;
}

// for vector printing
template<size_t N>
inline std::ostream &operator<<(std::ostream &out, const VecN<N> &a) {
    out << "[ ";
    for(size_t i=0; i < a.size();++i) {
        if(i>0) out << ", ";
        out << a[i];
    }
    out << " ]";
    return out;
}


} // close namespace pele
#endif
