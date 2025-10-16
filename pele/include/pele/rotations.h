#ifndef _pele_rotations_h_
#define _pele_rotations_h_

#include "pele/array.h"
#include "pele/vecn.h"

namespace pele{


/**
 * make a rotation matrix from an angle axis
 */
pele::MatrixNM<3,3> aa_to_rot_mat(pele::VecN<3> const & p);
pele::VecN<4> rot_mat_to_quaternion(pele::MatrixNM<3,3> const & mx);
pele::VecN<3> quaternion_to_aa(pele::VecN<4> const & q);
pele::VecN<4> aa_to_quaternion(pele::VecN<3> const & p);

inline pele::VecN<3> rot_mat_to_aa(pele::MatrixNM<3,3> const & mx)
{
    return pele::quaternion_to_aa(rot_mat_to_quaternion(mx));
}
pele::VecN<4> quaternion_multiply(pele::VecN<4> const & q1, pele::VecN<4> const & q2);
/**
 * change angle axis rotation p1 by the rotation p2
 */
inline pele::VecN<3> rotate_aa(pele::VecN<3> const & p1, pele::VecN<3> const & p2)
{
    return pele::quaternion_to_aa(
            pele::quaternion_multiply(
                pele::aa_to_quaternion(p2),
                pele::aa_to_quaternion(p1)
            )
    );
}

/**
 * compute the rotation matrix and it's derivatives from an angle axis vector if the rotation angle is very small
 */
void rot_mat_derivatives_small_theta(
        pele::VecN<3> const & p,
        pele::MatrixNM<3,3> & rmat,
        pele::MatrixNM<3,3> & drm1,
        pele::MatrixNM<3,3> & drm2,
        pele::MatrixNM<3,3> & drm3,
        bool with_grad);

/**
 * make a rotation matrix and it's derivatives from an angle axis
 */
void rot_mat_derivatives(
        pele::VecN<3> const & p,
        pele::MatrixNM<3,3> & rmat,
        pele::MatrixNM<3,3> & drm1,
        pele::MatrixNM<3,3> & drm2,
        pele::MatrixNM<3,3> & drm3);

} // close namespace pele
#endif
