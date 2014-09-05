/**
 * This is a partial c++ implementation of the rigid body potential and coordinate system.
 * This primarily just implements the functions necessary for potential calls.
 *
 *   - convert from rigid body coords to atomistic coords
 *   - convert an atomistic gradient into a gradient in the rb coordinate system.
 *
 * The next step is to compute distances between rigid body structures.  To do this we
 * will need to define a function which aligns the atomistic coordinates for fixed
 * rigid body centers of mass.  This will require
 *  - rot_mat_to_aa
 *  - each site must have a list of symmetries
 *  - rotate_aa()
 *
 *  Also need to replace the python zeroEv, which is composed of rotate, and align_path
 *  - align path is easy
 *  - rotate simply applies a rotation to a set of angle axis coords.  This will need to have access
 *    to the coords_adaptor
 *
 */
#ifndef _PELE_AATOPOLOGY_H_
#define _PELE_AATOPOLOGY_H_

#include <string>
#include <list>
#include <cmath>
#include <stdexcept>

#include "pele/array.h"
#include "pele/rotations.h"
#include "pele/base_potential.h"
#include "pele/vecn.h"
#include "pele/lowest_eig_potential.h"

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
    /**
     * the first dimension of the matrix
     *
     * note: if we make this const we can only use the assignment operator on
     * matrices with the same first dimension
     */
    size_t _dim2;
    HackyMatrix(size_t dim1, size_t dim2, dtype val=0)
        : pele::Array<dtype>(dim1 * dim2, val),
          _dim2(dim2)
    {}

    /**
     * wrap a pele::Array
     */
    HackyMatrix(pele::Array<double> v, size_t dim2)
        : pele::Array<dtype>(v),
          _dim2(dim2)
    {
        if (v.size() % dim2 != 0) {
            throw std::invalid_argument("v.size() is not divisible by dim2");
        }
    }

    /**
     * wrap an existing block of memory
     */
    HackyMatrix(double * data, size_t dim1, size_t dim2)
        : pele::Array<dtype>(data, dim1*dim2),
          _dim2(dim2)
    {}

    inline dtype const & operator()(size_t i, size_t j) const
    {
        return this->operator[](i * _dim2 + j);
    }
    inline dtype & operator()(size_t i, size_t j)
    {
        return this->operator[](i * _dim2 + j);
    }

    inline std::pair<size_t, size_t> shape() const
    {
        return std::pair<size_t, size_t>(this->size() / _dim2, _dim2);
    }

};

/**
 * multiply two matrices
 */
template<class dtype>
HackyMatrix<dtype> hacky_mat_mul(HackyMatrix<dtype> const & A, HackyMatrix<dtype> const & B)
{
    assert(A.shape().second == B.shape().first);
    size_t const L = A.shape().second;
    size_t const N = A.shape().first;
    size_t const M = B.shape().second;

    HackyMatrix<dtype> C(N, M, 0);
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

///**
// * multiply a matrix times an vector
// */
//template<class dtype>
//pele::Array<dtype> hacky_mat_mul(HackyMatrix<dtype> const & A, pele::Array<dtype> const & v)
//{
//    assert(A.shape().second == v.size());
//    size_t const L = A.shape().second;
//    size_t const n = A.shape().first;
//
//    pele::Array<dtype> C(n, 0);
//    for (size_t i = 0; i<n; ++i){
//        dtype val = 0;
//        for (size_t k = 0; k<L; ++k){
//            val += A(i,k) * v[k];
//        }
//        C(i) = val;
//    }
//    return C;
//}



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
        if (natoms != 0)
            throw std::runtime_error("coords adaptor doesn't support free atoms yet");
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

    pele::Array<double> get_rb_rotation(size_t isite)
    {
        if (_nrigid == 0) {
            // return empty array
            return pele::Array<double>();
        }
        if (isite > _nrigid) {
            throw std::invalid_argument("isite must be less than nrigid");
        }
        size_t const istart = 3*_nrigid + 3*isite;
        return _coords.view(istart, istart+3);
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
        if (_atom_positions.size() == 0 ) {
            throw std::invalid_argument("the atom positions must not have zero size");
        }
        if (_atom_positions.size() != _natoms * _ndim ) {
            throw std::invalid_argument("the length of atom_positions must be divisible by 3");
        }
    }

    inline size_t natoms() const { return _natoms; }

    /**
     * convert a center of mass and a angle axis rotation to a set of atomistic coordinates
     */
    pele::Array<double> to_atomistic(pele::Array<double> const com,
            pele::VecN<3> const & p)
    {
        assert(com.size() == _ndim);
        assert(p.size() == 3);
        auto rmat = pele::aa_to_rot_mat(p);
        Array<double> pos(_atom_positions.size());
        HackyMatrix<double> mpos(pos, _ndim);

        // in python this is:
        //      return com + np.dot(R, np.transpose(self.atom_positions)).transpose()
        for (size_t atom = 0; atom<_natoms; ++atom) {
            for (size_t j = 0; j<_ndim; ++j) {
                double val = com[j];
                for (size_t k = 0; k<_ndim; ++k) {
                    val += rmat(j,k) * _atom_positions_matrix(atom,k);
                }
                mpos(atom, j) = val;
            }
        }
        return pos;
    }

    /**
     * transform an atomistic gradient into a gradient in the
     * rigid body coordinates
     */
    void transform_grad(
            pele::VecN<3> const & p,
            pele::Array<double> const g,
            pele::VecN<3> & g_com,
            pele::VecN<3> & g_rot
            )
    {
        assert(g.size() == natoms() * 3);
        // view the array as a matrix
        HackyMatrix<double> gmat(g, 3);

        // compute the rotation matrix and derivatives
        pele::MatrixNM<3,3> rmat;
        pele::MatrixNM<3,3> drm1;
        pele::MatrixNM<3,3> drm2;
        pele::MatrixNM<3,3> drm3;
        rot_mat_derivatives(p, rmat, drm1, drm2, drm3);

        // do the center of mass coordinates
        for (size_t k=0; k<3; ++k) {
            double val = 0;
            for (size_t atom=0; atom < _natoms; ++atom) {
                val += gmat(atom,k);
            }
            g_com[k] = val;
        }

        // now do the rotations
        g_rot.assign(0);
        for (size_t atom=0; atom < _natoms; ++atom) {
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

    /**
     * transform an atomistic gradient into a gradient in the
     * rigid body coordinates
     *
     * This is simply a wrapper.  This copies the data into VecN objects,
     * calls transform_grad and copies it back.
     */
    void transform_grad(
            pele::Array<double> const & p,
            pele::Array<double> const g,
            pele::Array<double> & g_com,
            pele::Array<double> & g_rot
            )
    {
        pele::VecN<3> p_vec = p;
        pele::VecN<3> g_com_vec = g_com;
        pele::VecN<3> g_rot_vec = g_rot;
        transform_grad(p_vec, g, g_com_vec, g_rot_vec);
        /**
         * copy the data back into the arrays
         */
        for (size_t i = 0; i<3; ++i) {
            g_com[i] = g_com_vec[i];
            g_rot[i] = g_rot_vec[i];
        }

    }
};

///**
// * Angle axis topology
// *
// * An angle axis topology stores all topology information for an angle axis
// * system. The AATopology is composed of several angle axis sites,
// * which describe the shape of the angle axis site and each site carries a
// * position and orientation. Therefore, the length of the coordinate array
// * must be 6*number_of_sites.
// */
//class AATopology {
//};

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

    size_t number_of_non_rigid_atoms() { return 0; }

    CoordsAdaptor get_coords_adaptor(pele::Array<double> x)
    {
        return CoordsAdaptor(nrigid(), number_of_non_rigid_atoms(), x);
    }

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
        // view the atomistic coords as a matrix
        HackyMatrix<double> atomistic_mat(atomistic, 3);
        size_t istart = 0;
        for (size_t isite=0; isite<nrigid; ++isite) {
            VecN<3> psite = rb_rot.view(isite*3, isite*3+3);
            auto site_atom_positions = _sites[isite].to_atomistic(
                    rb_pos.view(isite*3, isite*3+3),
                    psite
                    );
            Array<double> atomistic_view(atomistic.view(istart, istart + site_atom_positions.size()));
            atomistic_view.assign(site_atom_positions);

            istart += site_atom_positions.size();
        }
        assert(istart == natoms_total() * 3);
        return atomistic;
    }

    /**
     * convert atomistic gradient into gradient in rigid body coordinates
     */
    void transform_gradient(pele::Array<double> rbcoords,
            pele::Array<double> grad, pele::Array<double> rbgrad)
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

    pele::VecN<3> align_angle_axis_vectors(pele::VecN<3> const & p1,
            pele::VecN<3> const & p2in)
    {
        pele::VecN<3> p2 = p2in;
        pele::VecN<3> n2, p2n;
        if (norm<3>(p2) < 1e-6) {
            if (norm<3>(p1) < 1e-6) {
                return p2;
            }
            n2 = p1;
            n2 *= 2. * M_PI / norm<3>(p1);
        } else {
            n2 = p2;
            n2 *= 2. * M_PI / norm<3>(p2);
        }

        while (true) {
            p2n = p2;
            p2n += n2;
            if (norm<3>(p2n - p1) > norm<3>(p2 - p1)) {
                break;
            }
            p2 = p2n;
        }

        while (true) {
            p2n = p2;
            p2n -= n2;
            if (norm<3>(p2n - p1) > norm<3>(p2 - p1)) {
                break;
            }
            p2 = p2n;
        }
        return p2;
    }

    /**
     * ensure a series of images are aligned with each other
     *
     * this simply aligns the angle axis vectors
     */
    void align_path(std::list<pele::Array<double> > path)
    {
        auto iter1 = path.begin();
        auto iter2 = path.begin();
        iter2++;
        while (iter2 != path.end()) {
            auto c1 = get_coords_adaptor(*iter1);
            auto c2 = get_coords_adaptor(*iter2);
            for (size_t isite = 0; isite < nrigid(); ++isite) {
                VecN<3> p1 = c1.get_rb_rotation(isite);
                pele::Array<double> p2 = c2.get_rb_rotation(isite);
                auto p2new = align_angle_axis_vectors(p1, p2);
                std::copy(p2new.begin(), p2new.end(), p2.begin());
            }
            ++iter1;
            ++iter2;
        }

    }

    /**
     * return a list of zero modes
     *
     * i.e. vectors corresponding to directions with zero curvature.
     * (these are not necessarily orthogonal)
     */
//    void get_zero_modes(pele::Array<double> const x,
//            std::vector<pele::Array<double> > & zev)
//    {
//        auto ca = get_coords_adaptor(x);
//        pele::Array<double> v(x.size(), 0);
//        auto cv = get_coords_adaptor(v);
//
//        // get the zero eigenvectors corresponding to translation
//        std::vector<pele::Array<double> > zev_t;
//        pele::zero_modes_translational(zev_t, nrigid(), 3);
//
//        for (auto const & v : zev_t) {
//            cv.get_rb_positions().assign(v);
//            zev.push_back(cv.get_coords().copy());
//        }
//
//        // get the zero eigenvectors corresponding to rotation
//        TransformAACluster transform(*this);
//        d = 1e-5
//        dx = x.copy()
//        transform.rotate(dx, rotations.aa2mx(np.array([d, 0, 0])))
//        self.align_path([x, dx])
//        dx -= x
//        dx /= np.linalg.norm(dx)
//
//        dy = x.copy()
//        transform.rotate(dy, rotations.aa2mx(np.array([0, d, 0])))
//        self.align_path([x, dy])
//        dy -= x
//        dy /= np.linalg.norm(dy)
//
//        dz = x.copy()
//        transform.rotate(dz, rotations.aa2mx(np.array([0, 0, d])))
//        self.align_path([x, dz])
//        dz -= x
//        dz /= np.linalg.norm(dz)
//
//        #print "Zero eigenvectors", zev
//        return zev + [dx, dy, dz]
//
//    }

};


class TransformPolicy {
public:
//    void translate(self, X, d) {
//        ''' translate the coordinates '''
//    }
    virtual ~TransformPolicy() {}

    /**
     *  apply rotation matrix mx for a rotation around the origin
     */
    virtual void rotate(pele::Array<double> x, pele::MatrixNM<3,3> const & mx) = 0;

//    def can_invert(self):
//        ''' returns True or False if an inversion can be performed'''
//
//    def invert(self, X):
//        ''' perform an inversion at the origin '''
//
//    def permute(self, X, perm):
//        ''' returns the permuted coordinates '''

};

class TransformAACluster : public TransformPolicy {
public:
    pele::RBTopology & m_topology;
    TransformAACluster(pele::RBTopology & topology)
        : m_topology(topology)
    {
    }
    virtual ~TransformAACluster() {}

    /**
     * apply a rotation to a set of rigid body coordinates
     */
    void rotate(pele::Array<double> x, pele::MatrixNM<3,3> const & mx)
    {
        auto ca = m_topology.get_coords_adaptor(x);
        if(m_topology.nrigid() > 0) {
            // rotate the center of mass positions by mx
            pele::HackyMatrix<double> rb_pos(ca.get_rb_positions(), 3);
            // make a HackyMatrix view of the transposed rotation matrix
            auto mxT = pele::transpose(mx);
            pele::HackyMatrix<double> mxT_view(mxT.data(), 3, 3);
            assert(mxT_view(0,1) == mx(1,0));
            // do the multiplication
            auto result = hacky_mat_mul(rb_pos, mxT_view);
            // copy the results back into the coordinates array
            std::cout << "result " << result << std::endl;
            rb_pos.assign(result);

            // rotate each aa rotation by mx
            VecN<3> dp = pele::rot_mat_to_aa(mx);
            auto rb_rot = ca.get_rb_rotations();
            for (size_t isite = 0; isite < m_topology.nrigid(); ++isite) {
                pele::Array<double> pview = rb_rot.view(isite*3, isite*3+3);
                VecN<3> p = pele::rotate_aa(pview, dp);
                // copy the vector back into pview
                std::copy(p.begin(), p.end(), pview.begin());
            }
        }
        if (m_topology.number_of_non_rigid_atoms() > 0) {
            throw std::runtime_error("non-rigid atoms is not yet supported");
//            ca.posAtom[:] = np.dot(mx, ca.posAtom.transpose()).transpose()
        }
    }

};


/**
 * potential wrapper for rigid body systems
 *
 * this converts rigid body coords to atomistic coords and passes the atomistic coords
 * to the potential function.  It also converts the atomistic gradient into a gradient
 * in the rigid body coordinate system.
 */
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

    inline double get_energy_gradient(pele::Array<double> rbcoords,
            pele::Array<double> rbgrad)
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
