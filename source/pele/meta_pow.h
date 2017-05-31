#ifndef _PELE_META_POW_H
#define _PELE_META_POW_H

/*
 *Usage: For integer (negative integer powers), one obtains x ^ N.
 *For half-integer (negative half-integer powers), one obtains sqrt(x ^ N).
 *Example: pos_int_pow<2>(x) is x ^ 2; neg_half_int_pow<5>(x) is x ^ (-2.5).
 *References:
 *Used here: general reference on template meta-programming and recursive template functions:
 *http://www.itp.phys.ethz.ch/education/hs12/programming_techniques
 *Also used here: reference on template meta-programming power function:
 *http://stackoverflow.com/questions/16443682/c-power-of-integer-template-meta-programming
 **/

#include <cmath> //for sqrt in half powers
#include <type_traits> //for static asserts

namespace pele{

template<class T, int N>
struct meta_pow{
    static T f(const T x)
    {
        return meta_pow<T, N - 1>::f(x) * x;
    }
};

template<class T>
struct meta_pow<T, 12>{
    static T f(const T x)
    {
        T x4 = meta_pow<T, 4>::f(x);
        return x4 * x4 * x4;
    }
};

template<class T>
struct meta_pow<T, 10>{
    static T f(const T x)
    {
        T x5 = meta_pow<T, 5>::f(x);
        return x5 * x5;
    }
};

template<class T>
struct meta_pow<T, 9>{
    static T f(const T x)
    {
        T x3 = x * x * x;
        return x3 * x3 * x3;
    }
};

template<class T>
struct meta_pow<T, 8>{
    static T f(const T x)
    {
        T x4 = meta_pow<T, 4>::f(x);
        return x4 * x4;
    }
};

template<class T>
struct meta_pow<T, 6>{
    static T f(const T x)
    {
        T x3 = x * x * x;
        return x3 * x3;
    }
};

template<class T>
struct meta_pow<T, 4>{
    static T f(const T x)
    {
        T x2 = x * x;
        return x2 * x2;
    }
};

template<class T>
struct meta_pow<T, 0>{
    static T f(const T x)
    {
        return T(1);
    }
};

/**
 * pow(x, N) where N >= 0, integer
 * usage: pos_int_pow<N>(x)
 */
template<int N, class T>
inline T pos_int_pow(const T x)
{
    static_assert(N >= 0, "illegal exponent input");
    return meta_pow<T, N>::f(x);
}

/**
 * pow(x, - N) where N >= 0, integer
 * usage: neg_int_pow<-N>(x)
 */
template<int N, class T>
inline T neg_int_pow(const T x)
{
    static_assert(N <= 0, "illegal exponent input");
    return T(1) / meta_pow<T, -N>::f(x);
}

/**
 * pow(x, N / 2) where N >= 0, integer
 * usage: pos_half_int_pow<N>(x)
 */
template<int N, class T>
inline T pos_half_int_pow(const T x)
{
    static_assert(N >= 0, "illegal exponent input");
    return std::sqrt(meta_pow<T, N>::f(x));
}

/**
 * pow(x, - N / 2) where N >= 0, integer
 * usage: neg_half_int_pow<-N>(x)
 */
template<int N, class T>
inline T neg_half_int_pow(const T x)
{
    static_assert(N <= 0, "illegal exponent input");
    return T(1) / std::sqrt(meta_pow<T, -N>::f(x));
}

}//namespace pele

#endif //#ifndef _PELE_META_POW_H
