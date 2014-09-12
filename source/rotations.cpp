#include "pele/rotations.h"
#include "pele/vecn.h"

namespace pele{


/**
 * make a rotation matrix from an angle axis
 */
MatrixNM<3,3> aa_to_rot_mat(pele::VecN<3> const & p)
{

    double theta2 = pele::dot<3>(p,p);
    MatrixNM<3,3> rm;
    if (theta2 < 1e-12) {
        MatrixNM<3,3> temp;
        rot_mat_derivatives_small_theta(p, rm, temp, temp, temp, false);
        return rm;
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
    pele::VecN<3> pn = p;
    pn*= theta;
    MatrixNM<3,3> e(0);
    e(0,1)  = -pn[2];
    e(0,2)  =  pn[1];
    e(1,2)  = -pn[0];
    e(1,0)  = -e(0,1);
    e(2,0)  = -e(0,2);
    e(2,1)  = -e(1,2);
    auto esq     = pele::dot<3,3,3>(e, e);

    // RM is calculated from Rodrigues' rotation formula [equation [1]
    // in the paper]
    // rm  = np.eye(3) + [1. - ct] * esq + st * e
    rm.assign(0);
    for (size_t i = 0; i < 3; ++i) rm(i,i) = 1; // identiy matrix
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            rm(i,j) += (1.-ct) * esq(i,j) + st * e(i,j);
        }
    }
    return rm;
}

pele::VecN<4> rot_mat_to_quaternion(pele::MatrixNM<3,3> const & mx)
{
    VecN<4> q(0);
    auto m = pele::transpose<3,3>(mx);
    double trace = m.trace();

    double s;
    if (trace > 0.){
        s = std::sqrt(trace+1.0) * 2.0;
        q[0] = 0.25 * s;
        q[1] = (m(1,2) - m(2,1)) / s;
        q[2] = (m(2,0) - m(0,2)) / s;
        q[3] = (m(0,1) - m(1,0)) / s;
    } else if ((m(0,0) > m(1,1)) and (m(0,0) > m(2,2))) {
        s = std::sqrt(1.0 + m(0,0) - m(1,1) - m(2,2)) * 2.0;
        q[0] = (m(1,2) - m(2,1)) / s;
        q[1] = 0.25 * s;
        q[2] = (m(1,0) + m(0,1)) / s;
        q[3] = (m(2,0) + m(0,2)) / s;
    } else if (m(1,1) > m(2,2)) {
        s = std::sqrt(1.0 + m(1,1) - m(0,0) - m(2,2)) * 2.0;
        q[0] = (m(2,0) - m(0,2)) / s;
        q[1] = (m(1,0) + m(0,1)) / s;
        q[2] = 0.25 * s;
        q[3] = (m(2,1) + m(1,2)) / s;
    } else {
        s = std::sqrt(1.0 + m(2,2) - m(0,0) - m(1,1)) * 2.0;
        q[0] = (m(0,1) - m(1,0)) / s;
        q[1] = (m(2,0) + m(0,2)) / s;
        q[2] = (m(2,1) + m(1,2)) / s;
        q[3] = 0.25 * s;
    }

    if (q[0] < 0) {
        q *= -1;
    }
    return q;
}

pele::VecN<3> quaternion_to_aa(pele::VecN<4> const & qin)
{
    double eps = 1e-6;
    VecN<3> p;
    VecN<4> q = qin;
    if (q[0] < 0.) {
        q *= -1;
    }
    if (q[0] > 1.0) {
        q /= std::sqrt(pele::dot<4>(q,q));
    }
    double theta = 2. * std::acos(q[0]);
    double s = std::sqrt(1. - q[0]*q[0]);
    for (size_t i = 0; i < 3; ++i) {
        p[i] = q[i+1];
    }
    if (s < eps) {
        p *= 2.;
    } else {
        p *= theta / s;
//        p = q[1:4] / s * theta
    }
    return p;
}

pele::VecN<4> aa_to_quaternion(pele::VecN<3> const & aa)
{
    static const double rot_epsilon = 1e-6;
    VecN<4> q;

    double thetah = 0.5 * pele::norm<3>(aa);
    q[0]  = std::cos(thetah);

    // do linear expansion for small epsilon
    VecN<3> v = aa;
    if (thetah < rot_epsilon) {
        v *= 0.5;
    } else {
        v *= 0.5 * std::sin(thetah) / thetah;
    }
    for (size_t i = 0; i < 3; ++i) {
        q[i+1] = v[i];
    }

    // make sure to have normal form
    if (q[0] < 0.0) q *= -1;
    return q;

}


pele::VecN<4> quaternion_multiply(pele::VecN<4> const & q0, pele::VecN<4> const & q1)
{
    VecN<4> q3;
    q3[0] = q0[0]*q1[0]-q0[1]*q1[1]-q0[2]*q1[2]-q0[3]*q1[3];
    q3[1] = q0[0]*q1[1]+q0[1]*q1[0]+q0[2]*q1[3]-q0[3]*q1[2];
    q3[2] = q0[0]*q1[2]-q0[1]*q1[3]+q0[2]*q1[0]+q0[3]*q1[1];
    q3[3] = q0[0]*q1[3]+q0[1]*q1[2]-q0[2]*q1[1]+q0[3]*q1[0];
    return q3;
}


/**
 * compute the rotation matrix and it's derivatives from an angle axis vector if the rotation angle is very small
 */
void rot_mat_derivatives_small_theta(
        pele::VecN<3> const & p,
        MatrixNM<3,3> & rmat,
        MatrixNM<3,3> & drm1,
        MatrixNM<3,3> & drm2,
        MatrixNM<3,3> & drm3,
        bool with_grad)
{
    double theta2 = pele::dot<3>(p, p);
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
 * make a rotation matrix and it's derivatives from an angle axis
 */
void rot_mat_derivatives(
        pele::VecN<3> const & p,
        MatrixNM<3,3> & rmat,
        MatrixNM<3,3> & drm1,
        MatrixNM<3,3> & drm2,
        MatrixNM<3,3> & drm3)
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

    double theta2 = pele::dot<3>(p,p);
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
    pele::VecN<3> pn = p;
    pn *= theta;
    MatrixNM<3,3> e(0);
    e(0,1)  = -pn[2];
    e(0,2)  =  pn[1];
    e(1,2)  = -pn[0];
    e(1,0)  = -e(0,1);
    e(2,0)  = -e(0,2);
    e(2,1)  = -e(1,2);
    auto esq     = pele::dot<3,3,3>(e, e);

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

    MatrixNM<3,3> de1(0.);
    de1(0,1) = p[2]*p[0]*theta3;
    de1(0,2) = -p[1]*p[0]*theta3;
    de1(1,2) = -(theta - p[0]*p[0]*theta3);
    de1(1,0) = -de1(0,1);
    de1(2,0) = -de1(0,2);
    de1(2,1) = -de1(1,2);

    MatrixNM<3,3> de2(0.);
    de2(0,1) = p[2]*p[1]*theta3;
    de2(0,2) = theta - p[1]*p[1]*theta3;
    de2(1,2) = p[0]*p[1]*theta3;
    de2(1,0) = -de2(0,1);
    de2(2,0) = -de2(0,2);
    de2(2,1) = -de2(1,2);

    MatrixNM<3,3> de3(0.);
    de3(0,1) = -(theta - p[2]*p[2]*theta3);
    de3(0,2) = -p[1]*p[2]*theta3;
    de3(1,2) = p[0]*p[2]*theta3;
    de3(1,0) = -de3(0,1);
    de3(2,0) = -de3(0,2);
    de3(2,1) = -de3(1,2);

    // compute the x derivative of the rotation matrix
    // drm1 = (st*pn[0]*esq + (1.-ct)*(de1.dot(e) + e.dot(de1))
    //         + ct*pn[0]*e + st*de1)
    auto de_e = pele::dot<3,3,3>(de1, e);
    auto ede_ = pele::dot<3,3,3>(e, de1);
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
    de_e = pele::dot<3,3,3>(de2, e);
    ede_ = pele::dot<3,3,3>(e, de2);
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
    de_e = pele::dot<3,3,3>(de3, e);
    ede_ = pele::dot<3,3,3>(e, de3);
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

} // namespace pele
