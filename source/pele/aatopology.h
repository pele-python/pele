#ifndef _PELE_AATOPOLOGY_H_
#define _PELE_AATOPOLOGY_H_

#include <string>
//#include <pair>
#include <cmath>

#include "pele/array.h"

namespace pele{

/**
 * this is a truly hacky implementation of a matrix.  please don't use it unless
 * you're being very careful
 */
template<class dtype>
class HackyMatrix : public pele::Array<dtype> {
public:
    size_t const _dim1;
    HackyMatrix(size_t dim1, size_t dim2, dtype val)
        : pele::Array<dtype>(dim1 * dim2, val),
          _dim1(dim1)
    {}
//    void set_dim1(size_t dim1)
//    {
//        if (this->size() % dim1 != 0) {
//            throw std::invalid_argument("size() is not divisible by dim1");
//        }
//        _dim1 = dim1;
//    }

    inline dtype const operator()(size_t i, size_t j) const
    {
        return this->operator[](i * _dim1 + j);
    }
    inline dtype & operator()(size_t i, size_t j)
    {
        return this->operator[](i * _dim1 + j);
    }

    inline std::pair<size_t, size_t> shape() const
    {
        return std::pair<size_t, size_t>(_dim1, this->size() / _dim1);
    }

};

template<class dtype>
HackyMatrix<dtype> hacky_mat_mul(HackyMatrix<dtype> const & A, HackyMatrix<dtype> const & B)
{
    assert(A.shape().second == B.shape().first);
    size_t const L = A.shape().second;
    size_t const n = A.shape().first;
    size_t const m = B.shape().second;

    HackyMatrix<dtype> C(n, m, 0);
    for (size_t i = 0; i<n; ++i){
        for (size_t j = 0; j<m; ++j){
            dtype val = 0;
            for (size_t k = 0; k<L; ++k){
                val += A[i,k] * B[k,j];
            }
            C[i,j] = val;
        }
    }
    return C;
}


/**
 * make a rotation matrix from an angle axis
 */
pele::Array<double> aa_to_rot_mat(pele::Array<double> const p)
{

    double theta2 = pele::norm(p);
//    if theta2 < 1e-12:
//        return _rot_mat_derivative_small_theta(p, with_grad)
    // Execute for the general case, where THETA dos not equal zero
    // Find values of THETA, CT, ST and THETA3
    double theta   = std::sqrt(theta2);
    double ct      = std::cos(theta);
    double st      = std::sin(theta);
    double theta3  = 1. / (theta2 * theta);

    // Set THETA to 1/THETA purely for convenience
    theta   = 1./theta;

    // Normalise p and construct the skew-symmetric matrix E
    // ESQ is calculated as the square of E
    auto pn = p.copy();
    pn*= theta;
    HackyMatrix<double> e(3, 3, 0);
    e(0,1)  = -pn[2];
    e(0,2)  =  pn[1];
    e(1,2)  = -pn[0];
    e(1,0)  = -e(0,1);
    e(2,0)  = -e(0,2);
    e(2,1)  = -e(1,2);
    auto esq     = hacky_mat_mul(e, e);

    // RM is calculated from Rodrigues' rotation formula [equation [1]
    // in the paper]
    // rm  = np.eye(3) + [1. - ct] * esq + st * e
    HackyMatrix<double> rm(3, 3, 0);
    for (size_t i = 0; i < 3; ++i) rm(i,i) = 1; // identiy matrix
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            rm(i,j) += (1.-ct) * esq(i,j) + st * e(i,j);
        }
    }
}

//class RBSite {

//};

/**
 * provide easy access to the different parts of a coordinates array
 *
 * the coords array will be filled as follows
 *
 * 0        -- 3*nrigid            : the center of mass of the rigid bodies
 * 3*nrigid -- 6*nrigid            : the rotations of the rigid bodies in angle axis coords
 * 6*nrigid -- 6*nrigid + 3*natoms : the positions of the non-rigid atoms (point masses)
 * ...      -- end                 : the last nlattice spaces are for the lattice degrees of freedom
 */
class CoordsAdaptor {
    size_t _nrigid;
    size_t _natoms;
    pele::Array<double> _coords;
    static const size_t _ndim = 3;

public:
    CoordsAdaptor(size_t nrigid, size_t natoms, pele::Array<double> coords)
        : _nrigid(nrigid),
          _natoms(natoms),
          _coords(coords)
    {
        if (coords.size() != (6*_nrigid + 3*_natoms)) {
            throw std::invalid_argument(std::string("coords has the wrong size ") + std::to_string(coords.size()));
        }
    }

    pele::Array<double> get_rb_positions()
    {
        if (_nrigid == 0) {
            // return empty array
            return pele::Array<double>();
        }
        return _coords.view(0, 3*_nrigid);
    }

    pele::Array<double> get_rb_rotations()
    {
        if (_nrigid == 0) {
            // return empty array
            return pele::Array<double>();
        }
        return _coords.view(3*_nrigid, 6*_nrigid);
    }

    pele::Array<double> get_atom_positions()
    {
        if (_natoms == 0) {
            // return empty array
            return pele::Array<double>();
        }

        return _coords.view(6*_nrigid, 6*_nrigid + 3*_natoms);
    }

};

}

#endif
