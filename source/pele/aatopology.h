/**
 * This is a partial c++ implementation of the rigid body potential and coordinate system.
 * This primarily just implements the functions necessary for potential calls.
 *
 *   - convert from rigid body coords to atomistic coords
 *   - convert an atomistic gradient into a gradient in the rb coordinate system.
 *
 */
#ifndef _PELE_AATOPOLOGY_H_
#define _PELE_AATOPOLOGY_H_

#include <string>
//#include <pair>
#include <cmath>
#include <stdexcept>

#include "pele/array.h"
#include "pele/base_potential.h"

namespace pele{

/**
 * this is a truly hacky implementation of a matrix.  please don't use it unless
 * you're being very careful
 *
 * it is a simply wrapper for pele::Array, so a pele array can be wrapped to
 * act as a matrix temporarily.  The idea is to redo somthing like the reshape() function
 * in numpy.
 */
template<class dtype>
class HackyMatrix : public pele::Array<dtype> {
public:
    size_t _dim1;
    HackyMatrix(size_t dim1, size_t dim2, dtype val=0)
        : pele::Array<dtype>(dim1 * dim2, val),
          _dim1(dim1)
    {}

    /**
     * wrap a pele::Array
     */
    HackyMatrix(pele::Array<double> v, size_t dim1)
        : pele::Array<dtype>(v),
          _dim1(dim1)
    {
        if (v.size() % dim1 != 0) {
            throw std::invalid_argument("v.size() is not divisible by dim1");
        }
    }

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
                val += A(i,k) * B(k,j);
            }
            C(i,j) = val;
        }
    }
    return C;
}

template<class dtype>
pele::Array<dtype> hacky_mat_mul(HackyMatrix<dtype> const & A, pele::Array<dtype> const & v)
{
    assert(A.shape().second == v.size());
    size_t const L = A.shape().second;
    size_t const n = A.shape().first;

    pele::Array<dtype> C(n, 0);
    for (size_t i = 0; i<n; ++i){
        dtype val = 0;
        for (size_t k = 0; k<L; ++k){
            val += A(i,k) * v[k];
        }
        C(i) = val;
    }
    return C;
}

/**
 * compute the rotation matrix and it's derivatives from an angle axis vector if the rotation angle is very small
 */
void rot_mat_derivatives_small_theta(
        pele::Array<double> const p,
        HackyMatrix<double> rmat,
        HackyMatrix<double> drm1,
        HackyMatrix<double> drm2,
        HackyMatrix<double> drm3,
        bool with_grad)
{
    double theta2 = dot(p, p);
    if (theta2 > 1e-2) {
        throw std::invalid_argument("theta must be small");
    }
    // Execute if the angle of rotation is zero
    // In this case the rotation matrix is the identity matrix
    rmat.assign(0);
    for (size_t i = 0; i<3; ++i) rmat(i,i) = 1; // identity matrix


    // vr274> first order corrections to rotation matrix
    rmat(0,1) = -p[2];
    rmat(1,0) = p[2];
    rmat(0,2) = p[1];
    rmat(2,0) = -p[1];
    rmat(1,2) = -p[0];
    rmat(2,1) = p[0];

    // If derivatives do not need to found, we're finished
    if (not with_grad) {
        return;
    }

    // hk286 - now up to the linear order in theta
    drm1.assign(0);
    drm1(0,0)    = 0.0;
    drm1(0,1)    = p[1];
    drm1(0,2)    = p[2];
    drm1(1,0)    = p[1];
    drm1(1,1)    = -2.0*p[0];
    drm1(1,2)    = -2.0;
    drm1(2,0)    = p[2];
    drm1(2,1)    = 2.0;
    drm1(2,2)    = -2.0*p[0];
    drm1 *= 0.5;

    drm2.assign(0);
    drm2(0,0)    = -2.0*p[1];
    drm2(0,1)    = p[0];
    drm2(0,2)    = 2.0;
    drm2(1,0)    = p[0];
    drm2(1,1)    = 0.0;
    drm2(1,2)    = p[2];
    drm2(2,0)    = -2.0;
    drm2(2,1)    = p[2];
    drm2(2,2)    = -2.0*p[1];
    drm2 *= 0.5;

    drm3.assign(0);
    drm3(0,0)    = -2.0*p[2];
    drm3(0,1)    = -2.0;
    drm3(0,2)    = p[0];
    drm3(1,0)    = 2.0;
    drm3(1,1)    = -2.0*p[2];
    drm3(1,2)    = p[1];
    drm3(2,0)    = p[0];
    drm3(2,1)    = p[1];
    drm3(2,2)    = 0.0;
    drm3 *= 0.5;
}


/**
 * make a rotation matrix from an angle axis
 */
pele::HackyMatrix<double> aa_to_rot_mat(pele::Array<double> const p)
{

    double theta2 = pele::dot(p,p);
    if (theta2 < 1e-12) {
        pele::HackyMatrix<double> rmat(3,3);
        pele::HackyMatrix<double> temp(3,3);
        rot_mat_derivatives_small_theta(p, rmat, temp, temp, temp, false);
        return rmat;
    }
    // Execute for the general case, where THETA dos not equal zero
    // Find values of THETA, CT, ST and THETA3
    double theta   = std::sqrt(theta2);
    double ct      = std::cos(theta);
    double st      = std::sin(theta);

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
    return rm;
}

/**
 * make a rotation matrix and it's derivatives from an angle axis
 */
void rot_mat_derivatives(
        pele::Array<double> const p,
        HackyMatrix<double> rmat,
        HackyMatrix<double> drm1,
        HackyMatrix<double> drm2,
        HackyMatrix<double> drm3)
{
    assert(p.size() == 3);
    if (rmat.shape() != std::pair<size_t, size_t>(3,3)) {
        throw std::invalid_argument("rmat matrix has the wrong size");
    }
    if (drm1.shape() != std::pair<size_t, size_t>(3,3)) {
        throw std::invalid_argument("drm1 matrix has the wrong size");
    }
    if (drm2.shape() != std::pair<size_t, size_t>(3,3)) {
        throw std::invalid_argument("drm2 matrix has the wrong size");
    }
    if (drm3.shape() != std::pair<size_t, size_t>(3,3)) {
        throw std::invalid_argument("drm2 matrix has the wrong size");
    }

    double theta2 = pele::dot(p,p);
    if (theta2 < 1e-12) {
        return rot_mat_derivatives_small_theta(p, rmat, drm1, drm2, drm3, true);
    }
    // Execute for the general case, where THETA dos not equal zero
    // Find values of THETA, CT, ST and THETA3
    double theta   = std::sqrt(theta2);
    double ct      = std::cos(theta);
    double st      = std::sin(theta);
    double theta3  = 1./(theta2*theta);


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

    // compute the rotation matrix
    // RM is calculated from Rodrigues' rotation formula [equation [1]
    // in the paper]
    // rm  = np.eye(3) + [1. - ct] * esq + st * e
    rmat.assign(0.);
    for (size_t i = 0; i < 3; ++i) rmat(i,i) = 1; // identiy matrix
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            rmat(i,j) += (1.-ct) * esq(i,j) + st * e(i,j);
        }
    }

    HackyMatrix<double> de1(3,3,0.);
    de1(0,1) = p[2]*p[0]*theta3;
    de1(0,2) = -p[1]*p[0]*theta3;
    de1(1,2) = -(theta - p[0]*p[0]*theta3);
    de1(1,0) = -de1(0,1);
    de1(2,0) = -de1(0,2);
    de1(2,1) = -de1(1,2);

    HackyMatrix<double> de2(3,3,0.);
    de2(0,1) = p[2]*p[1]*theta3;
    de2(0,2) = theta - p[1]*p[1]*theta3;
    de2(1,2) = p[0]*p[1]*theta3;
    de2(1,0) = -de2(0,1);
    de2(2,0) = -de2(0,2);
    de2(2,1) = -de2(1,2);

    HackyMatrix<double> de3(3,3,0.);
    de3(0,1) = -(theta - p[2]*p[2]*theta3);
    de3(0,2) = -p[1]*p[2]*theta3;
    de3(1,2) = p[0]*p[2]*theta3;
    de3(1,0) = -de3(0,1);
    de3(2,0) = -de3(0,2);
    de3(2,1) = -de3(1,2);

    // compute the x derivative of the rotation matrix
    // drm1 = (st*pn[0]*esq + (1.-ct)*(de1.dot(e) + e.dot(de1))
    //         + ct*pn[0]*e + st*de1)
    auto de_e = hacky_mat_mul(de1, e);
    auto ede_ = hacky_mat_mul(e, de1);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            drm1(i,j) = (st * pn[0] * esq(i,j)
                         + (1.-ct) * (de_e(i,j) + ede_(i,j))
                         + ct * pn[0] * e(i,j)
                         + st * de1(i,j)
                         );
        }
    }

    // compute the y derivative of the rotation matrix
    de_e = hacky_mat_mul(de2, e);
    ede_ = hacky_mat_mul(e, de2);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            drm2(i,j) = (st * pn[1] * esq(i,j)
                         + (1.-ct) * (de_e(i,j) + ede_(i,j))
                         + ct * pn[1] * e(i,j)
                         + st * de2(i,j)
                         );
        }
    }

    // compute the z derivative of the rotation matrix
    de_e = hacky_mat_mul(de3, e);
    ede_ = hacky_mat_mul(e, de3);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            drm3(i,j) = (st * pn[2] * esq(i,j)
                         + (1.-ct) * (de_e(i,j) + ede_(i,j))
                         + ct * pn[2] * e(i,j)
                         + st * de3(i,j)
                         );
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
    size_t _nrigid;  /** the number of rigid bodies */
    size_t _natoms;  /** the number of non-rigid point particles */
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

    pele::Array<double> get_coords() { return _coords; }

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

/**
 * represent a single rigid body
 */
class RigidFragment {
    static const size_t _ndim = 3;
    pele::Array<double> _atom_positions;
    pele::HackyMatrix<double> _atom_positions_matrix;
    size_t _natoms;

public:
    RigidFragment(pele::Array<double> atom_positions)
    : _atom_positions(atom_positions.copy()),
      _atom_positions_matrix(_atom_positions, _ndim),
      _natoms(_atom_positions.size() / _ndim)
    {
//        std::cout << "atoms positions " << _atom_positions << "\n";
        if (_atom_positions.size() == 0 ) {
            throw std::invalid_argument("the atoms positions must not have zero size");
        }
        if (_atom_positions.size() != _natoms * _ndim ) {
            throw std::invalid_argument("the length of atom_positions must be divisible by 3");
        }
    }

    size_t natoms() const { return _natoms; }

    /**
     * convert a center of mass and a angle axis rotation to a set of atomistic coordinates
     */
    pele::Array<double> to_atomistic(pele::Array<double> const com, pele::Array<double> const p)
    {
        assert(com.size() == _ndim);
        assert(p.size() == 3);
        auto rmat = pele::aa_to_rot_mat(p);
        Array<double> pos(_atom_positions.size());
        HackyMatrix<double> mpos(pos, _ndim);

        // return com + np.dot(R, np.transpose(self.atom_positions)).transpose()
        for (size_t atom = 0; atom<_natoms; ++atom) {
            for (size_t j = 0; j<_ndim; ++j) {
                double val = com[j];
                for (size_t k = 0; k<_ndim; ++k) {
                    val += rmat(j,k) * _atom_positions_matrix(atom,k);
                }
                mpos(atom, j) = val;
            }
        }
//        std::cout << _atom_positions << "\n";
//        std::cout << pos << "\n";
//        std::cout << rmat << "\n";
//        std::cout << mpos << "\n";
        return pos;
    }

    /**
     * transform an atomistic gradient into a gradient in the
     * rigid body coordinates
     */
    void transform_grad(
            pele::Array<double> const p,
            pele::Array<double> const g,
            pele::Array<double> g_com,
            pele::Array<double> g_rot
            )
    {
        assert(p.size() == 3);
        assert(g.size() == natoms() * 3);
        assert(g_com.size() == 3);
        assert(g_rot.size() == 3);
        HackyMatrix<double> gmat(g, 3);

        // compute the rotation matrix derivatives
        HackyMatrix<double> rmat(3,3);
        HackyMatrix<double> drm1(3,3);
        HackyMatrix<double> drm2(3,3);
        HackyMatrix<double> drm3(3,3);
//        std::cout << rmat.shape().first << " " << rmat.shape().second << "\n";
        rot_mat_derivatives(p, rmat, drm1, drm2, drm3);

        // do the center of mass coordinates
        for (size_t k=0; k<3; ++k) {
            double val = 0;
            for (size_t atom=0; atom<natoms(); ++atom) {
                val += gmat(atom,k);
            }
            g_com[k] = val;
        }

        // now do the rotations
        g_rot.assign(0);
        for (size_t atom=0; atom<natoms(); ++atom) {
            double val1 = 0;
            double val2 = 0;
            double val3 = 0;
            for (size_t i=0; i<3; ++i) {
                for (size_t j=0; j<3; ++j) {
                    val1 += gmat(atom,i) * drm1(i,j) * _atom_positions_matrix(atom,j);
                    val2 += gmat(atom,i) * drm2(i,j) * _atom_positions_matrix(atom,j);
                    val3 += gmat(atom,i) * drm3(i,j) * _atom_positions_matrix(atom,j);
                }
            }
            g_rot[0] += val1;
            g_rot[1] += val2;
            g_rot[2] += val3;
        }

    }
};

/**
 * represent a collection of rigid bodies
 */
class RBTopology {
    std::vector<RigidFragment> _sites;
    size_t _natoms_total;
//    bool _finalized;
//    std::vector<std::vector<double> > _atom_indices;

public:
    RBTopology()
        : _natoms_total(0)//, _finalized(false)
    {}

    void add_site(Array<double> atom_positions)
    {
        _sites.push_back(RigidFragment(atom_positions));
    }

    void finalize()
    {
        _natoms_total = 0;
        for (auto & rf : _sites) {
//            _atom_indices.push_back(std::vector<double>(rf.natoms()));
//            for (size_t i = 0; i<rf.natoms(); ++i) {
//                _atom_indices.back()[i] = _natoms + i;
//            }
            _natoms_total += rf.natoms();
        }
    }

    /**
     * number of rigid bodies
     */
    size_t nrigid() const { return _sites.size(); }
    /**
     * return the total number of atoms in the atomistic representation
     */
    size_t natoms_total() const { return _natoms_total; }

    Array<double> to_atomistic(Array<double> rbcoords)
    {
        if (natoms_total() == 0) {
            finalize();
        }
        if ( rbcoords.size() != nrigid() * 6 ) {
            throw std::invalid_argument("rbcoords has the wrong size");
        }

        size_t const nrigid = _sites.size();
        CoordsAdaptor ca(nrigid, 0, rbcoords);
        auto rb_pos = ca.get_rb_positions();
        auto rb_rot = ca.get_rb_rotations();
        Array<double> atomistic(3 * natoms_total());
        HackyMatrix<double> atomistic_mat(atomistic, 3);
        size_t istart = 0;
        for (size_t isite=0; isite<nrigid; ++isite) {
            auto site_atom_positions = _sites[isite].to_atomistic(
                    rb_pos.view(isite*3, isite*3+3),
                    rb_rot.view(isite*3, isite*3+3)
                    );
            Array<double> atomistic_view(atomistic.view(istart, istart + site_atom_positions.size()));
            atomistic_view.assign(site_atom_positions);

//            HackyMatrix<double> site_atom_positions_mat(site_atom_positions, 3);
//            for (size_t iatom = 0; iatom<_sites[isite].natoms(); ++iatom) {
//                for (size_t k = 0; k<3; ++k) {
//                    size_t const atom = _atom_indices[isite][iatom];
//                    atomistic_mat(atom, k) = site_atom_positions_mat(iatom, k)
//                }
//            }

            istart += site_atom_positions.size();
        }
        assert(istart == natoms_total() * 3);
        return atomistic;
    }

    /**
     * convert atomistic gradient into gradient in rigid body coordinates
     */
    void transform_gradient(pele::Array<double> rbcoords, pele::Array<double> grad, pele::Array<double> rbgrad)
    {
        if (natoms_total() == 0) {
            finalize();
        }
        if ( rbcoords.size() != nrigid() * 6 ) {
            throw std::invalid_argument("rbcoords has the wrong size");
        }
        if (grad.size() != natoms_total() * 3) {
            throw std::invalid_argument("grad has the wrong size");
        }
        if (rbgrad.size() != rbcoords.size()) {
            throw std::invalid_argument("rbgrad has the wrong size");
        }

        CoordsAdaptor ca(nrigid(), 0, rbcoords);
        pele::Array<double> coords_rot(ca.get_rb_rotations());
//        pele::Array<double> rbgrad(rbcoords.size());
        CoordsAdaptor rbgrad_ca(nrigid(), 0, rbgrad);
        HackyMatrix<double> g_com(rbgrad_ca.get_rb_positions(), 3);
        HackyMatrix<double> g_rot(rbgrad_ca.get_rb_rotations(), 3);

        size_t istart = 0;
        for (size_t isite=0; isite<nrigid(); ++isite) {
            size_t const site_ndof = _sites[isite].natoms() * 3;
//            std::cout << grad.size() << " " << istart << " " << site_ndof << " " << istart + site_ndof << "\n";
            Array<double> g_site     = grad.view      (istart, istart + site_ndof);
            Array<double> p          = coords_rot.view(isite*3, isite*3 + 3);
            Array<double> g_com_site = g_com.view     (isite*3, isite*3 + 3);
            Array<double> g_rot_site = g_rot.view     (isite*3, isite*3 + 3);
            _sites[isite].transform_grad(p, g_site, g_com_site, g_rot_site);
            istart += site_ndof;
        }
    }
};

class RBPotentialWrapper : public BasePotential {
    std::shared_ptr<BasePotential> potential_;
    RBTopology topology_;
public:

    RBPotentialWrapper(std::shared_ptr<BasePotential> potential)
        : potential_(potential)
    {}

    inline void add_site(pele::Array<double> atom_positions)
    {
        topology_.add_site(atom_positions);
    }

    inline double get_energy(pele::Array<double> rbcoords)
    {
        auto x = topology_.to_atomistic(rbcoords);
        return potential_->get_energy(x);
    }

    inline double get_energy_gradient(pele::Array<double> rbcoords, pele::Array<double> rbgrad)
    {
        auto x = topology_.to_atomistic(rbcoords);
        pele::Array<double> grad_atomistic(topology_.natoms_total() * 3);
        double e = potential_->get_energy_gradient(x, grad_atomistic);
        topology_.transform_gradient(rbcoords, grad_atomistic, rbgrad);
        return e;
    }


};

}

#endif
