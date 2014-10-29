/**
 * This is a partial c++ implementation of the tools needed to interact with
 * systems of rigid bodies.  This is not a complete reimplementation, only the
 * parts that were too slow in python were implemented here.
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
 * This is a very minimal implementation of a matrix.  It's primary function is
 * to act as a wrapper for pele::Array, so a pele array can be act as a matrix
 * temporarily.  The idea is to redo somthing like the reshape() function in
 * numpy.
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
 * 6*nrigid -- 6*nrigid + 3*natoms : the positions of the non-rigid atoms (point masses) (not yet supported)
 * ...      -- end                 : the last nlattice spaces are for the lattice degrees of freedom (not yet supported)
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

    /**
     * return the full coords array
     */
    pele::Array<double> get_coords() { return _coords; }

    /**
     * return a view of the rigid body centers of mass
     */
    pele::Array<double> get_rb_positions()
    {
        if (_nrigid == 0) {
            // return empty array
            return pele::Array<double>();
        }
        return _coords.view(0, 3*_nrigid);
    }

    /**
     * return a view of the center of mass coords of a specific rigid body
     */
    pele::Array<double> get_rb_position(size_t isite)
    {
        if (_nrigid == 0) {
            // return empty array
            return pele::Array<double>();
        }
        if (isite > _nrigid) {
            throw std::invalid_argument("isite must be less than nrigid");
        }
        size_t const istart = 3*isite;
        return _coords.view(istart, istart+3);
    }

    /**
     * return a view of the rigid angle axis rotations
     */
    pele::Array<double> get_rb_rotations()
    {
        if (_nrigid == 0) {
            // return empty array
            return pele::Array<double>();
        }
        return _coords.view(3*_nrigid, 6*_nrigid);
    }

    /**
     * return a view of the angle axis rotation of a specific rigid body
     */
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

// forward definition of RBTopology needed for TrasnformAACluster
class RBTopology;

/**
 * This is the base class from which all Transform Policies should be derived
 */
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

/**
 * Routines to apply transformations to a rigid body cluster
 */
class TransformAACluster : public TransformPolicy {
public:
    pele::RBTopology * m_topology;
    TransformAACluster(pele::RBTopology * topology)
        : m_topology(topology)
    { }
    virtual ~TransformAACluster() {}

    /**
     * apply a rotation to a set of rigid body coordinates
     */
    void rotate(pele::Array<double> x, pele::MatrixNM<3,3> const & mx);
//    inline void rotate(pele::Array<double> x, pele::Array<double> mx)
//    {
//        return rotate(x, pele::MatrixNM<3,3>(mx));
//    }
};

class MeasureAngleAxisCluster {
public:
    pele::RBTopology * m_topology;
    MeasureAngleAxisCluster(pele::RBTopology * topology)
        : m_topology(topology)
    { }

    /**
     * align the rotations so that the atomistic coordinates will be in best alignment
     */
    void align(pele::Array<double> const x1, pele::Array<double> x2);


};


/**
 * represent a single rigid body
 */
class RigidFragment {
    static const size_t _ndim = 3;
    pele::Array<double> _atom_positions;
    pele::HackyMatrix<double> _atom_positions_matrix;
    size_t _natoms;

    double m_M; // total mass of the angle axis site
    double m_W; // sum of all weights
    pele::VecN<3> m_cog; // center of gravity
    pele::MatrixNM<3,3> m_S; // weighted tensor of gyration S_ij = \sum m_i x_i x_j
    pele::MatrixNM<3,3> m_inversion; // matrix that applies the appropriate inversion
    bool m_can_invert;

    // a list of rotations that leave the rigid body unchanged.
    std::vector<pele::MatrixNM<3,3> > m_symmetry_rotations;


public:
    RigidFragment(pele::Array<double> atom_positions,
            Array<double> cog,
            double M,
            double W,
            Array<double> S,
            Array<double> inversion, bool can_invert)
    : _atom_positions(atom_positions.copy()),
      _atom_positions_matrix(_atom_positions, _ndim),
      _natoms(_atom_positions.size() / _ndim),
      m_M(M),
      m_W(W),
      m_cog(cog),
      m_S(S),
      m_inversion(inversion),
      m_can_invert(can_invert)
    {
        if (_atom_positions.size() == 0 ) {
            throw std::invalid_argument("the atom positions must not have zero size");
        }
        if (_atom_positions.size() != _natoms * _ndim ) {
            throw std::invalid_argument("the length of atom_positions must be divisible by 3");
        }
    }

    /**
     * return the number of atoms in the rigid body
     */
    inline size_t natoms() const { return _natoms; }

    /**
     * add a symmetry rotation
     */
    inline void add_symmetry_rotation(pele::Array<double> R)
    {
        m_symmetry_rotations.push_back(R);
    }

    /**
     * access the vector of symmetry rotations
     */
    inline std::vector<pele::MatrixNM<3,3> > const & get_symmetry_rotations() const
    {
        return m_symmetry_rotations;
    }

    /**
     * convert a center of mass and a angle axis rotation to a set of atomistic coordinates
     */
    pele::Array<double> to_atomistic(pele::Array<double> const com,
            pele::VecN<3> const & p);

    /**
     * transform an atomistic gradient into a gradient in the
     * rigid body coordinates
     */
    void transform_grad(
            pele::VecN<3> const & p,
            pele::Array<double> const g,
            pele::VecN<3> & g_com,
            pele::VecN<3> & g_rot
            );

    /**
     * transform an atomistic gradient into a gradient in the
     * rigid body coordinates
     *
     * This is simply a wrapper.  This copies the data into VecN objects,
     * calls transform_grad and copies it back.
     */
    void transform_grad(
            pele::Array<double> const p,
            pele::Array<double> const g,
            pele::Array<double> g_com,
            pele::Array<double> g_rot
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

    /**
     * return the shortest vector from com1 to com2
     *
     * this could be replaced by periodic distances for instance
     */
    inline pele::VecN<3> get_smallest_rij(pele::VecN<3> const & com1, pele::VecN<3> const & com2) const
    {
        return com2 - com1;
    }


    // sn402 version: overloaded to include a boxlength vector
    inline pele::VecN<3> get_smallest_rij_bulk(pele::VecN<3> const & com1, pele::VecN<3> const & com2, pele::VecN<3> const & boxvec) const
    {
        pele::VecN<3> dx = com2 - com1;
        pele::VecN<3> iboxvec;

        for (int i=0;i<3;i++)
        {
        	iboxvec[i]=1/boxvec[i];
        	dx[i] -= boxvec[i] * round(dx[i]*iboxvec[i]);
        }
        return dx;
    }

    /**
     * compute the squared distance between two configurations of the rigid fragment
     */
    double distance_squared(pele::VecN<3> const & com1, pele::VecN<3> const & p1,
            pele::VecN<3> const & com2, pele::VecN<3> const & p2) const;

    // sn402 overloaded header
    double distance_squared_bulk(pele::VecN<3> const & com1, pele::VecN<3> const & p1,
            pele::VecN<3> const & com2, pele::VecN<3> const & p2, pele::VecN<3> const & boxvec) const;




    void distance_squared_grad(pele::VecN<3> const & com1, pele::VecN<3> const & p1,
            pele::VecN<3> const & com2, pele::VecN<3> const & p2,
            VecN<3> & g_M, VecN<3> & g_P
            ) const;

    // sn402 overloaded header
    void distance_squared_grad_bulk(pele::VecN<3> const & com1, pele::VecN<3> const & p1,
            pele::VecN<3> const & com2, pele::VecN<3> const & p2,
            VecN<3> & g_M, VecN<3> & g_P, pele::VecN<3> const & boxvec
            ) const;


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

public:
    RBTopology()
        : _natoms_total(0)
    {}

    void add_site(RigidFragment const & site)
    {
        _sites.push_back(site);
        _natoms_total += site.natoms();
    }

    /**
     * provide access to the vector of rigid fragments
     */
    std::vector<RigidFragment> const & get_sites() const { return _sites; };

    /**
     * return the number of rigid bodies
     */
    size_t nrigid() const { return _sites.size(); }

    /**
     * return the total number of atoms in the atomistic representation
     */
    size_t natoms_total() const { return _natoms_total; }

    size_t number_of_non_rigid_atoms() const { return 0; }

    /**
     * return an already constructed CoordsAdaptor object
     */
    CoordsAdaptor get_coords_adaptor(pele::Array<double> x) const
    {
        return CoordsAdaptor(nrigid(), number_of_non_rigid_atoms(), x);
    }

    /**
     * convert rigid body coordinates to atomistic coordinates
     */
    Array<double> to_atomistic(Array<double> rbcoords);

    /**
     * convert atomistic gradient into gradient in rigid body coordinates
     */
    void transform_gradient(pele::Array<double> rbcoords,
            pele::Array<double> grad, pele::Array<double> rbgrad);

    /**
     * align two angle axis vectors
     *
     * perform symmetry operations on p2 to minimize the distance with p1
     */
    pele::VecN<3> align_angle_axis_vectors(pele::VecN<3> const & p1,
            pele::VecN<3> const & p2in);

    /**
     * Ensure the angle axis rotations of two structures are aligned with each other.
     *
     * x1 will remain unchanged, only modify x2
     */
    void align_all_angle_axis_vectors(pele::Array<double> x1,
            pele::Array<double> x2);

    /**
     * ensure a series of images are aligned with each other
     *
     * this simply aligns the angle axis vectors
     */
    void align_path(std::list<pele::Array<double> > path);

    /**
     * return a list of zero modes
     *
     * i.e. vectors corresponding to directions with zero curvature.
     * (these are not necessarily orthogonal)
     */
    void get_zero_modes(pele::Array<double> const x,
            std::vector<pele::Array<double> > & zev)
    {
        auto ca = get_coords_adaptor(x);
        pele::Array<double> v(x.size(), 0);
        auto cv = get_coords_adaptor(v);

        // get the zero eigenvectors corresponding to translation
        std::vector<pele::Array<double> > zev_t;
        pele::zero_modes_translational(zev_t, nrigid(), 3);

        for (auto const & v : zev_t) {
            cv.get_rb_positions().assign(v);
            zev.push_back(cv.get_coords().copy());
        }

        // get the zero eigenvectors corresponding to rotation
        TransformAACluster transform(this);
        double d = 1e-5;
        pele::VecN<3> v3;
        pele::Array<double> delta(x.size());

        // do rotations around the x y and z axes
        for (size_t i = 0; i < 3; ++i) {
            delta.assign(x);
            v3.assign(0);
            v3[i] = d;
            transform.rotate(delta, pele::aa_to_rot_mat(v3));
            align_all_angle_axis_vectors(x, delta);
            delta -= x;
            delta /= norm(delta);
            zev.push_back(delta.copy());
        }
    }

    /**
     * return the squared distance between two configurations
     */
    double distance_squared(pele::Array<double> const x1, pele::Array<double> const x2) const
    {
        double d_sq = 0;
        auto ca1 = get_coords_adaptor(x1);
        auto ca2 = get_coords_adaptor(x2);
        for (size_t isite = 0; isite < nrigid(); ++isite) {
            d_sq += _sites[isite].distance_squared(
                    ca1.get_rb_position(isite),
                    ca1.get_rb_rotation(isite),
                    ca2.get_rb_position(isite),
                    ca2.get_rb_rotation(isite)
                );
        }
        return d_sq;
    }

    // sn402 version: overloaded to include a boxlength vector.
    double distance_squared_bulk(pele::Array<double> const x1, pele::Array<double> const x2, pele::VecN<3> const boxvec) const
    {
        double d_sq = 0;
        auto ca1 = get_coords_adaptor(x1);
        auto ca2 = get_coords_adaptor(x2);
        for (size_t isite = 0; isite < nrigid(); ++isite) {
            d_sq += _sites[isite].distance_squared_bulk(
                    ca1.get_rb_position(isite),
                    ca1.get_rb_rotation(isite),
                    ca2.get_rb_position(isite),
                    ca2.get_rb_rotation(isite),
                    boxvec
                );
        }
        return d_sq;
    }

    /**
     * Calculate gradient with respect to x1 for the squared distance
     *
     * used to compute the spring force on x1 to x2
     */
    void distance_squared_grad(pele::Array<double> const x1, pele::Array<double> const x2,
            pele::Array<double> grad
            ) const
    {
        if (grad.size() != x1.size()) {
            throw std::runtime_error("grad has the wrong size");
        }
        grad.assign(0);
        auto ca1 = get_coords_adaptor(x1);
        auto ca2 = get_coords_adaptor(x2);
        auto ca_spring = get_coords_adaptor(grad);

        // first distance for sites only
        for (size_t isite=0; isite<nrigid(); ++isite) {
            pele::VecN<3> g_M, g_P;
            _sites[isite].distance_squared_grad(
                    ca1.get_rb_position(isite),
                    ca1.get_rb_rotation(isite),
                    ca2.get_rb_position(isite),
                    ca2.get_rb_rotation(isite),
                    g_M, g_P);
            auto spring_com = ca_spring.get_rb_position(isite);
            std::copy(g_M.begin(), g_M.end(), spring_com.begin());
            auto spring_rot = ca_spring.get_rb_rotation(isite);
            std::copy(g_P.begin(), g_P.end(), spring_rot.begin());
        }
    }



    // sn402 overloaded version to include boxvec
    void distance_squared_grad_bulk(pele::Array<double> const x1, pele::Array<double> const x2,
    						   pele::Array<double> grad, pele::VecN<3> const boxvec
        					  ) const
    {
    	if (grad.size() != x1.size())
    	{
    		throw std::runtime_error("grad has the wrong size");
    	}
    	grad.assign(0);
    	auto ca1 = get_coords_adaptor(x1);
    	auto ca2 = get_coords_adaptor(x2);
    	auto ca_spring = get_coords_adaptor(grad);

    	// first distance for sites only
    	for (size_t isite=0; isite<nrigid(); ++isite)
    	{
    		pele::VecN<3> g_M, g_P;
    		_sites[isite].distance_squared_grad_bulk(
                ca1.get_rb_position(isite),
                ca1.get_rb_rotation(isite),
                ca2.get_rb_position(isite),
                ca2.get_rb_rotation(isite),
                g_M, g_P, boxvec);
    		auto spring_com = ca_spring.get_rb_position(isite);
    		std::copy(g_M.begin(), g_M.end(), spring_com.begin());
    		auto spring_rot = ca_spring.get_rb_rotation(isite);
    		std::copy(g_P.begin(), g_P.end(), spring_rot.begin());
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
    std::shared_ptr<RBTopology> topology_;
public:

    RBPotentialWrapper(std::shared_ptr<BasePotential> potential,
            std::shared_ptr<RBTopology> top)
        : potential_(potential),
          topology_(top)

    {}

    inline double get_energy(pele::Array<double> rbcoords)
    {
        auto x = topology_->to_atomistic(rbcoords);
        return potential_->get_energy(x);
    }

    inline double get_energy_gradient(pele::Array<double> rbcoords,
            pele::Array<double> rbgrad)
    {
        auto x = topology_->to_atomistic(rbcoords);
        pele::Array<double> grad_atomistic(topology_->natoms_total() * 3);
        double e = potential_->get_energy_gradient(x, grad_atomistic);
        topology_->transform_gradient(rbcoords, grad_atomistic, rbgrad);
        return e;
    }
};

}

#endif
