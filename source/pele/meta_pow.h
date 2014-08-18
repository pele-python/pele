#ifndef _PELE_META_POW_H
#define _PELE_META_POW_H

namespace pele{

#include <cmath> //for sqrt in half powers

template<class T, int N>
struct meta_pow{
    static T f(const T x)
    {
        return meta_pow<T, N - 1>::f(x) * x;
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
    return meta_pow<T, N>::f(x);
}

/**
 * pow(x, - N) where N >= 0, integer
 * usage: neg_int_pow<-N>(x)
 */
template<int N, class T>
inline T neg_int_pow(const T x)
{
    return T(1) / meta_pow<T, -N>::f(x);
}

/**
 * pow(x, N / 2) where N >= 0, integer
 * usage: pos_half_int_pow<N>(x)
 */
template<int N, class T>
inline T pos_half_int_pow(const T x)
{
    return std::sqrt(meta_pow<T, N>::f(x));
}

/**
 * pow(x, - N / 2) where N >= 0, integer
 * usage: neg_half_int_pow<-N>(x)
 */
template<int N, class T>
inline T neg_half_int_pow(const T x)
{
    return T(1) / std::sqrt(meta_pow<T, -N>::f(x));
}

}//namespace pele

#endif //#ifndef _PELE_META_POW_H
